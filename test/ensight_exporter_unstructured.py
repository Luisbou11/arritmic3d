"""
ensight_exporter_unstructured.py
=================================
Exportador Arritmic3D -> EnSight Gold unstructured (HEXA8).

Fase 2: geometría (.geo) + variables estáticas y temporales (.scl) + .case
completo con secciones VARIABLE y TIME.

Uso
---
    from ensight_exporter_unstructured import EnsightGoldWriter

    writer = EnsightGoldWriter.from_pyvista_grid(
        case_dir         = "output/ensight_unstructured",
        base_name        = "ventricle",
        tissue           = tissue,
        grid             = grid,
        static_grid_vars = ["Cell_type"],
        static_data      = {"TissueType": np.array(v_type, dtype=float)},
        variables        = ["State", "APD", "CV"],
        binary           = True,   # False para ASCII
    )

    writer.write_timestep(tissue, time_ms=20.0)
    writer.finalize()

Lógica de la máscara de nodos
------------------------------
    1. Parámetro explícito node_mask (toma precedencia).
    2. static_data["TissueType"] != 0
    3. static_data["Cell_type"] != -1
    4. ValueError si no se puede inferir.

Indexado
--------
    global_id = z * nx * ny + y * nx + x   (z más lento, x más rápido)
    local_id (EnSight): 1-based.

Vértices HEXA8 (EnSight Gold / VTK_HEXAHEDRON)
-----------------------------------------------
    n1 = (i,   j,   k  )    n5 = (i,   j,   k+1)
    n2 = (i+1, j,   k  )    n6 = (i+1, j,   k+1)
    n3 = (i+1, j+1, k  )    n7 = (i+1, j+1, k+1)
    n4 = (i,   j+1, k  )    n8 = (i,   j+1, k+1)
"""

from __future__ import annotations

import io
import os
from typing import Sequence

import numpy as np

# ---------------------------------------------------------------------------
# Tabla de variables soportadas: nombre público -> método de CardiacTissue
# ---------------------------------------------------------------------------
_VARIABLE_MAP: dict[str, str] = {
    "State":  "GetStates",
    "APD":    "GetAPD",
    "CV":     "GetCV",
    "DI":     "GetDI",
    "LastDI": "GetLastDI",
    "LAT":    "GetLAT",
    "Beat":   "GetBeat",
    "AP":     "GetAP",
}

_ENSIGHT_MAX_DESC = 79
_ENSIGHT_STR_LEN  = 80


def _truncate(text: str, max_len: int = _ENSIGHT_MAX_DESC) -> str:
    return text[:max_len]


def _str80(text: str) -> bytes:
    """Codifica text como campo de exactamente 80 bytes relleno con espacios."""
    encoded = text.encode("ascii", errors="replace")
    return encoded[:_ENSIGHT_STR_LEN].ljust(_ENSIGHT_STR_LEN)


class EnsightGoldWriter:
    """
    Exportador Arritmic3D -> EnSight Gold unstructured (HEXA8).
    Escribe geometría + variables estáticas y temporales en ASCII.
    """

    # ------------------------------------------------------------------
    # Constructor principal
    # ------------------------------------------------------------------
    def __init__(
        self,
        case_dir: str,
        base_name: str,
        tissue,
        variables: Sequence[str] | None = None,
        static_data: dict[str, np.ndarray] | None = None,
        x_coords: list[float] | None = None,
        y_coords: list[float] | None = None,
        z_coords: list[float] | None = None,
        binary: bool = True,
        node_mask: np.ndarray | None = None,
    ) -> None:
        self.case_dir  = os.path.abspath(case_dir)
        self.base_name = base_name
        self.binary    = binary
        self.variables = list(variables) if variables is not None else []
        self._static_data: dict[str, np.ndarray] = static_data or {}
        self._timestep_index  = 0
        self._timestep_values: list[float] = []

        os.makedirs(self.case_dir, exist_ok=True)

        if x_coords is not None:
            self._x_coords = np.asarray(x_coords, dtype=float)
            self._y_coords = np.asarray(y_coords, dtype=float)
            self._z_coords = np.asarray(z_coords, dtype=float)
        else:
            self._ensure_coords(tissue)

        # Infer node mask (explicit > TissueType > Cell_type > error)
        if node_mask is not None:
            self._node_mask = np.asarray(node_mask, dtype=bool)
        elif "TissueType" in self._static_data:
            self._node_mask = self._static_data["TissueType"] != 0
        elif "Cell_type" in self._static_data:
            self._node_mask = self._static_data["Cell_type"] != -1
        else:
            raise ValueError(
                "No se puede inferir la máscara de nodos válidos. "
                "Proporciona 'TissueType' o 'Cell_type' en static_data, "
                "o pasa node_mask explícitamente al constructor."
            )

        self._build_geometry(tissue)
        self._write_geometry(tissue)

        # Write static variables and update .case with VARIABLE section
        for var_name, data in self._static_data.items():
            self._write_static_variable(var_name, data)

        self._update_case_file_variables()

        mode = "binario" if binary else "ASCII"
        print(
            f"[EnsightGoldWriter-Unstructured] Inicializado en '{self.case_dir}' (modo {mode}).",
            flush=True,
        )

    # ------------------------------------------------------------------
    # Constructor alternativo desde grid pyvista
    # ------------------------------------------------------------------
    @classmethod
    def from_pyvista_grid(
        cls,
        case_dir: str,
        base_name: str,
        tissue,
        grid,
        variables: Sequence[str] | None = None,
        static_grid_vars: Sequence[str] | None = None,
        static_data: dict[str, np.ndarray] | None = None,
        binary: bool = True,
        node_mask: np.ndarray | None = None,
    ) -> "EnsightGoldWriter":
        """
        Crea un EnsightGoldWriter extrayendo las coordenadas reales del grid pyvista.

        static_grid_vars: nombres de arrays en grid.point_data que se exportan
            como variables estáticas, p.ej. ["Cell_type"].
        static_data: arrays precalculados adicionales, p.ej. {"TissueType": ...}.
        En caso de nombre duplicado prevalece static_data.
        node_mask: máscara booleana explícita de nodos válidos (opcional).
        """
        x_coords = np.unique(grid.points[:, 0]).tolist()
        y_coords = np.unique(grid.points[:, 1]).tolist()
        z_coords = np.unique(grid.points[:, 2]).tolist()

        merged: dict[str, np.ndarray] = {}
        if static_grid_vars:
            for name in static_grid_vars:
                merged[name] = np.asarray(grid.point_data[name], dtype=float)
        if static_data:
            merged.update(static_data)

        return cls(
            case_dir   = case_dir,
            base_name  = base_name,
            tissue     = tissue,
            variables  = variables,
            static_data = merged,
            x_coords   = x_coords,
            y_coords   = y_coords,
            z_coords   = z_coords,
            binary     = binary,
            node_mask  = node_mask,
        )

    # ------------------------------------------------------------------
    # API pública: timesteps
    # ------------------------------------------------------------------
    def write_timestep(self, tissue, time_ms: float) -> None:
        """Escribe un snapshot de todas las variables temporales al tiempo time_ms."""
        self._timestep_index += 1
        self._timestep_values.append(float(time_ms))
        for var_name in self.variables:
            getter = getattr(tissue, _VARIABLE_MAP[var_name])
            data_full = np.asarray(getter(), dtype=float)
            data_filtered = data_full[self._node_mask]
            self._write_variable(var_name, data_filtered, self._timestep_index)

    def finalize(self) -> None:
        """Actualiza el .case con la sección TIME final e imprime resumen."""
        self._write_case_file()

        n_ts     = self._timestep_index
        n_vars   = len(self.variables)
        n_static = len(self._static_data)
        n_scl    = n_static + n_vars * n_ts
        vars_str = ", ".join(self.variables) if self.variables else "(ninguna)"

        print(
            f"[EnsightGoldWriter-Unstructured] Simulación completada.\n"
            f"Timesteps escritos: {n_ts}\n"
            f"Variables: {vars_str} ({n_vars} variables)\n"
            f"Directorio de salida: {self.case_dir}\n"
            f"Ficheros: {self.base_name}.geo, {self.base_name}.case, "
            f"{self.base_name}_*.scl ({n_scl} ficheros de variables)",
            flush=True,
        )

    # ------------------------------------------------------------------
    # Construcción de la geometría en memoria
    # ------------------------------------------------------------------
    def _build_geometry(self, tissue) -> None:
        nx = tissue.GetSizeX()
        ny = tissue.GetSizeY()
        nz = tissue.GetSizeZ()
        n_total = nx * ny * nz

        mask = self._node_mask
        if len(mask) != n_total:
            raise ValueError(
                f"node_mask tiene {len(mask)} elementos, esperado {n_total} "
                f"(nx={nx}, ny={ny}, nz={nz})."
            )

        # Coordinates for every global_id = z*nx*ny + y*nx + x.
        ZZ, YY, XX = np.meshgrid(
            np.arange(nz, dtype=np.intp),
            np.arange(ny, dtype=np.intp),
            np.arange(nx, dtype=np.intp),
            indexing="ij",
        )
        all_x = self._x_coords[XX.ravel()]
        all_y = self._y_coords[YY.ravel()]
        all_z = self._z_coords[ZZ.ravel()]

        # Keep only nodes inside the mask
        self._node_x = all_x[mask]
        self._node_y = all_y[mask]
        self._node_z = all_z[mask]
        n_nodes = int(mask.sum())

        # global_id -> local_id (1-indexed; 0 = invalid / outside mask)
        self._global_to_local = np.zeros(n_total, dtype=np.int32)
        self._global_to_local[mask] = np.arange(1, n_nodes + 1, dtype=np.int32)

        # HEXA8 connectivity: iterate over all potential cells
        ci = np.arange(nx - 1, dtype=np.intp)
        cj = np.arange(ny - 1, dtype=np.intp)
        ck = np.arange(nz - 1, dtype=np.intp)
        CI, CJ, CK = np.meshgrid(ci, cj, ck, indexing="ij")
        CI = CI.ravel()
        CJ = CJ.ravel()
        CK = CK.ravel()

        def _gid(ix, jy, kz):
            return kz * (nx * ny) + jy * nx + ix

        G = np.stack([
            _gid(CI,     CJ,     CK    ),  # n1
            _gid(CI + 1, CJ,     CK    ),  # n2
            _gid(CI + 1, CJ + 1, CK    ),  # n3
            _gid(CI,     CJ + 1, CK    ),  # n4
            _gid(CI,     CJ,     CK + 1),  # n5
            _gid(CI + 1, CJ,     CK + 1),  # n6
            _gid(CI + 1, CJ + 1, CK + 1),  # n7
            _gid(CI,     CJ + 1, CK + 1),  # n8
        ], axis=1)

        # Keep only cells whose 8 corners are all valid
        valid = np.all(mask[G], axis=1)
        G_valid = G[valid]

        self._hex_connectivity = self._global_to_local[G_valid]  # shape: (n_hex, 8)

        n_hex = len(self._hex_connectivity)
        ratio = n_nodes / n_total if n_total > 0 else 0.0
        print(
            f"[EnsightGoldWriter-Unstructured] "
            f"n_nodes_total={n_total:,}, n_nodes_kept={n_nodes:,}, "
            f"n_hex={n_hex:,}, ratio={ratio:.4f}",
            flush=True,
        )

    # ------------------------------------------------------------------
    # Geometría (.geo) — despachador
    # ------------------------------------------------------------------
    def _write_geometry(self, tissue) -> None:
        if self.binary:
            self._write_geometry_binary(tissue)
        else:
            self._write_geometry_ascii(tissue)

    def _write_geometry_binary(self, tissue) -> None:
        n_nodes = len(self._node_x)
        n_hex   = len(self._hex_connectivity)

        geo_path = os.path.join(self.case_dir, f"{self.base_name}.geo")

        with open(geo_path, "wb") as f:
            f.write(_str80("C Binary"))
            f.write(_str80("EnSight Gold Geometry File"))
            f.write(_str80("Arritmic3D cardiac tissue export (unstructured)"))
            f.write(_str80("node id off"))
            f.write(_str80("element id off"))

            f.write(_str80("part"))
            f.write(np.int32(1).tobytes())
            f.write(_str80("CardiacTissue"))

            f.write(_str80("coordinates"))
            f.write(np.int32(n_nodes).tobytes())

            f.write(np.asarray(self._node_x, dtype=np.float32).tobytes())
            f.write(np.asarray(self._node_y, dtype=np.float32).tobytes())
            f.write(np.asarray(self._node_z, dtype=np.float32).tobytes())

            f.write(_str80("hexa8"))
            f.write(np.int32(n_hex).tobytes())
            f.write(self._hex_connectivity.astype(np.int32).tobytes())

        print(
            f"[EnsightGoldWriter-Unstructured] Geometría escrita (binario): '{geo_path}'",
            flush=True,
        )

    def _write_geometry_ascii(self, tissue) -> None:
        n_nodes = len(self._node_x)
        n_hex   = len(self._hex_connectivity)

        geo_path = os.path.join(self.case_dir, f"{self.base_name}.geo")

        with open(geo_path, "w") as f:
            f.write(_truncate("EnSight Gold Geometry File") + "\n")
            f.write(_truncate("Arritmic3D cardiac tissue export (unstructured)") + "\n")
            f.write("node id off\n")
            f.write("element id off\n")

            f.write("part\n")
            f.write(f"{1:10d}\n")
            f.write(_truncate("CardiacTissue") + "\n")

            f.write("coordinates\n")
            f.write(f"{n_nodes:10d}\n")

            # Component-major: all X, then all Y, then all Z — one value per line
            self._write_float_block_np(f, self._node_x)
            self._write_float_block_np(f, self._node_y)
            self._write_float_block_np(f, self._node_z)

            # Connectivity: 8 local IDs per line, %10d
            f.write("hexa8\n")
            f.write(f"{n_hex:10d}\n")
            conn_fmt = " ".join(["%10d"] * 8)
            np.savetxt(f, self._hex_connectivity, fmt=conn_fmt)

        print(
            f"[EnsightGoldWriter-Unstructured] Geometría escrita (ASCII): '{geo_path}'",
            flush=True,
        )

    @staticmethod
    def _write_float_block_np(f, arr: np.ndarray) -> None:
        """Escribe un array en formato EnSight ASCII: %12.5e, 1 valor por línea."""
        buf = io.StringIO()
        np.savetxt(buf, arr, fmt="%12.5e")
        f.write(buf.getvalue())

    # ------------------------------------------------------------------
    # Variables estáticas (.scl sin timestep) — despachador + implementaciones
    # ------------------------------------------------------------------
    def _write_static_variable(self, var_name: str, data: np.ndarray) -> None:
        """Filtra data con _node_mask y escribe la variable estática."""
        data_filtered = np.asarray(data, dtype=float)[self._node_mask]
        if self.binary:
            self._write_static_variable_binary(var_name, data_filtered)
        else:
            self._write_static_variable_ascii(var_name, data_filtered)

    def _write_static_variable_binary(self, var_name: str, data: np.ndarray) -> None:
        filename = f"{self.base_name}_{var_name}.scl"
        filepath = os.path.join(self.case_dir, filename)
        with open(filepath, "wb") as f:
            f.write(_str80(f"Arritmic3D {var_name} (static)"))
            f.write(_str80("part"))
            f.write(np.int32(1).tobytes())
            f.write(_str80("coordinates"))
            f.write(np.asarray(data, dtype=np.float32).tobytes())
        print(
            f"[EnsightGoldWriter-Unstructured] Variable estática escrita (binario): '{filepath}'",
            flush=True,
        )

    def _write_static_variable_ascii(self, var_name: str, data: np.ndarray) -> None:
        """
        Escribe <base_name>_<var_name>.scl para una variable estática.
        data ya viene filtrado (longitud n_nodes).
        """
        filename = f"{self.base_name}_{var_name}.scl"
        filepath = os.path.join(self.case_dir, filename)
        with open(filepath, "w") as f:
            f.write(_truncate(f"Arritmic3D {var_name} (static)") + "\n")
            f.write("part\n")
            f.write(f"{1:10d}\n")
            f.write("coordinates\n")
            buf = io.StringIO()
            np.savetxt(buf, data, fmt="%12.5e")
            f.write(buf.getvalue())
        print(
            f"[EnsightGoldWriter-Unstructured] Variable estática escrita (ASCII): '{filepath}'",
            flush=True,
        )

    # ------------------------------------------------------------------
    # Variables temporales (.scl por timestep) — despachador + implementaciones
    # ------------------------------------------------------------------
    def _write_variable(self, var_name: str, data, timestep_idx: int) -> None:
        """data ya viene filtrado (longitud n_nodes)."""
        if self.binary:
            self._write_variable_binary(var_name, data, timestep_idx)
        else:
            self._write_variable_ascii(var_name, data, timestep_idx)

    def _write_variable_binary(self, var_name: str, data, timestep_idx: int) -> None:
        filename = f"{self.base_name}_{var_name}.{timestep_idx:05d}.scl"
        filepath = os.path.join(self.case_dir, filename)
        with open(filepath, "wb") as f:
            f.write(_str80(f"Arritmic3D {var_name}"))
            f.write(_str80("part"))
            f.write(np.int32(1).tobytes())
            f.write(_str80("coordinates"))
            f.write(np.asarray(data, dtype=np.float32).tobytes())

    def _write_variable_ascii(self, var_name: str, data, timestep_idx: int) -> None:
        """
        Escribe <base_name>_<var_name>.<NNNNN>.scl para el timestep dado.
        data ya viene filtrado (longitud n_nodes).
        """
        filename = f"{self.base_name}_{var_name}.{timestep_idx:05d}.scl"
        filepath = os.path.join(self.case_dir, filename)
        with open(filepath, "w") as f:
            f.write(_truncate(f"Arritmic3D {var_name}") + "\n")
            f.write("part\n")
            f.write(f"{1:10d}\n")
            f.write("coordinates\n")
            buf = io.StringIO()
            np.savetxt(buf, np.asarray(data, dtype=float), fmt="%12.5e")
            f.write(buf.getvalue())

    # ------------------------------------------------------------------
    # Fichero .case — con secciones VARIABLE y TIME
    # ------------------------------------------------------------------
    def _update_case_file_variables(self) -> None:
        """Actualiza el .case para incluir la sección VARIABLE (estáticas + temporales)."""
        self._write_case_file()

    def _write_case_file(self) -> None:
        """Escribe el fichero .case con FORMAT, GEOMETRY, VARIABLE y TIME."""
        case_path = os.path.join(self.case_dir, f"{self.base_name}.case")

        with open(case_path, "w") as f:
            f.write("FORMAT\n")
            f.write("type: ensight gold\n")
            f.write("\n")
            f.write("GEOMETRY\n")
            f.write(f"model: {self.base_name}.geo\n")

            if self._static_data or self.variables:
                f.write("\n")
                f.write("VARIABLE\n")
                # Variables estáticas: referencia directa al .scl
                for var_name in self._static_data:
                    filename = f"{self.base_name}_{var_name}.scl"
                    desc = _truncate(var_name, 49)
                    f.write(f"scalar per node:     {desc:<24s} {filename}\n")
                # Variables temporales: wildcard sobre timesteps, time set 1
                for var_name in self.variables:
                    wildcard = f"{self.base_name}_{var_name}.*****.scl"
                    desc = _truncate(var_name, 49)
                    f.write(f"scalar per node:  1  {desc:<24s} {wildcard}\n")

            if self._timestep_values:
                f.write("\n")
                f.write("TIME\n")
                f.write("time set:              1\n")
                f.write(f"number of steps:       {len(self._timestep_values)}\n")
                f.write("filename start number: 1\n")
                f.write("filename increment:    1\n")
                f.write("filename number width: 5\n")
                f.write("time values:\n")
                for t in self._timestep_values:
                    f.write(f"{t:12.5e}\n")

    # ------------------------------------------------------------------
    # Utilidades internas
    # ------------------------------------------------------------------
    def _ensure_coords(self, tissue) -> None:
        """Coordenadas enteras de fallback cuando no se dispone de grid pyvista."""
        nx = tissue.GetSizeX()
        ny = tissue.GetSizeY()
        nz = tissue.GetSizeZ()
        self._x_coords = np.arange(nx, dtype=float)
        self._y_coords = np.arange(ny, dtype=float)
        self._z_coords = np.arange(nz, dtype=float)
