"""
ensight_exporter.py
===================
Módulo de exportación Arritmic3D -> EnSight Gold (ASCII).

Genera el fichero de geometría (.geo), los ficheros de variables escalares
(.scl) y el descriptor (.case) usando el formato 'block rectilinear', que es
el más eficiente para mallas rectilineares y está soportado por ParaView y
otras herramientas.

Uso
---
    from ensight_exporter import EnsightGoldWriter

    # Con espaciado explícito
    writer = EnsightGoldWriter(
        case_dir  = "/ruta/salida",
        base_name = "cardiac",
        tissue    = tissue,
        variables = ["State", "APD", "LAT"],
        x_coords  = [...],
        y_coords  = [...],
        z_coords  = [...],
    )

    # Con grid pyvista
    writer = EnsightGoldWriter.from_pyvista_grid(
        case_dir  = "/ruta/salida",
        base_name = "cardiac",
        tissue    = tissue,
        variables = ["State", "APD", "LAT"],
        grid      = grid,
    )

    # En el bucle de simulación
    writer.write_timestep(tissue, time_ms=20.0)
    writer.write_timestep(tissue, time_ms=40.0)
    writer.finalize()

Estructura de ficheros generada
--------------------------------
    <case_dir>/
        <base_name>.case
        <base_name>.geo
        <base_name>_State_00001.scl
        <base_name>_APD_00001.scl
        ...

Variables disponibles
---------------------
    "State", "APD", "CV", "DI", "LastDI", "LAT", "Beat", "AP"
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


def _truncate(text: str, max_len: int = _ENSIGHT_MAX_DESC) -> str:
    return text[:max_len]


class EnsightGoldWriter:
    """
    Exportador Arritmic3D -> EnSight Gold ASCII.
    Genera .case, .geo (block rectilinear) y .scl (escalares por timestep).
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
    ) -> None:
        self.case_dir   = os.path.abspath(case_dir)
        self.base_name  = base_name
        self.variables  = list(variables) if variables is not None else []
        self._static_data: dict[str, np.ndarray] = static_data or {}
        self._timestep_index  = 0
        self._timestep_values: list[float] = []

        os.makedirs(self.case_dir, exist_ok=True)

        if x_coords is not None:
            self._x_coords = list(x_coords)
            self._y_coords = list(y_coords)
            self._z_coords = list(z_coords)
        else:
            self._ensure_coords(tissue)

        self._write_geometry(tissue)
        for var_name, data in self._static_data.items():
            self._write_static_variable(var_name, data)
        self._write_case_file()

        print(
            f"[EnsightGoldWriter] Inicializado en '{self.case_dir}'.",
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
    ) -> "EnsightGoldWriter":
        """
        Crea un EnsightGoldWriter extrayendo las coordenadas reales del grid pyvista.

        static_grid_vars: nombres de arrays en grid.point_data que se exportan
            como variables estáticas, p.ej. ["Cell_type"].
        static_data: arrays precalculados adicionales que no provienen del grid,
            p.ej. {"TissueType": np.array(v_type)}.
        Ambos parámetros se combinan; en caso de nombre duplicado prevalece static_data.
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
            case_dir    = case_dir,
            base_name   = base_name,
            tissue      = tissue,
            variables   = variables,
            static_data = merged,
            x_coords    = x_coords,
            y_coords    = y_coords,
            z_coords    = z_coords,
        )

    # ------------------------------------------------------------------
    # API pública: timesteps
    # ------------------------------------------------------------------
    def write_timestep(self, tissue, time_ms: float) -> None:
        """Escribe un snapshot de todas las variables al tiempo time_ms."""
        self._timestep_index += 1
        self._timestep_values.append(float(time_ms))
        for var_name in self.variables:
            getter = getattr(tissue, _VARIABLE_MAP[var_name])
            data   = getter()
            self._write_variable(var_name, data, self._timestep_index)
        self._write_case_file()

    def finalize(self) -> None:
        """Actualiza el .case con el número final de timesteps."""
        self._write_case_file()

    # ------------------------------------------------------------------
    # Geometría (.geo)
    # ------------------------------------------------------------------
    def _write_geometry(self, tissue) -> None:
        """
        Escribe <base_name>.geo en formato EnSight Gold ASCII.

        Usa 'block rectilinear': solo requiere las coordenadas únicas de
        cada eje (ni + nj + nk valores), sin expandir la malla completa.
        ParaView y otras herramientas lo soportan correctamente.

        Formato de coordenadas: %12.5e, 1 valor por línea.
        """
        nx = tissue.GetSizeX()
        ny = tissue.GetSizeY()
        nz = tissue.GetSizeZ()

        geo_path = os.path.join(self.case_dir, f"{self.base_name}.geo")

        with open(geo_path, "w") as f:
            # Cabeceras obligatorias
            f.write(_truncate("EnSight Gold Geometry File") + "\n")
            f.write(_truncate("Arritmic3D cardiac tissue export") + "\n")
            f.write("node id off\n")
            f.write("element id off\n")

            # Parte (ID base-1 obligatorio)
            f.write("part\n")
            f.write(f"{1:10d}\n")
            f.write(_truncate("CardiacTissue") + "\n")

            # block rectilinear: solo coordenadas únicas por eje
            f.write("block rectilinear\n")

            # Dimensiones ni nj nk en una sola línea
            f.write(f"{nx:10d}{ny:10d}{nz:10d}\n")

            # Coordenadas X (ni valores), Y (nj valores), Z (nk valores)
            self._write_float_block_np(f, np.asarray(self._x_coords))
            self._write_float_block_np(f, np.asarray(self._y_coords))
            self._write_float_block_np(f, np.asarray(self._z_coords))

        print(f"[EnsightGoldWriter] Geometría escrita: '{geo_path}'", flush=True)

    @staticmethod
    def _write_float_block_np(f, arr: np.ndarray) -> None:
        """
        Escribe un numpy array en formato EnSight: %12.5e, 1 valor por línea.

        BUG CORREGIDO: la versión anterior escribía 6 valores por línea
        (np.savetxt con reshape(-1, 6)). La especificación EnSight Gold exige
        exactamente 1 valor por línea para block rectilinear ('E12.5  1/line').
        Con 6/línea el lector tomaba solo el primer valor de cada línea, lo que
        desalineaba los tres ejes: los nodos X recibían valores de Y, Y recibía
        valores de Z, y Z quedaba a cero. El resultado era una malla desplazada
        y colapsada visible en ParaView.
        """
        buf = io.StringIO()
        np.savetxt(buf, arr, fmt="%12.5e")
        f.write(buf.getvalue())

    # ------------------------------------------------------------------
    # Variables escalares (.scl)
    # ------------------------------------------------------------------
    def _write_static_variable(self, var_name: str, data: np.ndarray) -> None:
        """
        Escribe una variable estática (sin timestep): <base_name>_<var_name>.scl.
        Se referencia en el .case sin time set, por lo que EnSight/ParaView la
        trata como constante en el tiempo.
        """
        filename = f"{self.base_name}_{var_name}.scl"
        filepath = os.path.join(self.case_dir, filename)
        with open(filepath, "w") as f:
            f.write(_truncate(f"Arritmic3D {var_name} (static)") + "\n")
            f.write("part\n")
            f.write(f"{1:10d}\n")
            f.write("block\n")
            buf = io.StringIO()
            np.savetxt(buf, np.asarray(data, dtype=float), fmt="%12.5e")
            f.write(buf.getvalue())
        print(f"[EnsightGoldWriter] Variable estatica escrita: '{filepath}'", flush=True)

    def _write_variable(self, var_name: str, data, timestep_idx: int) -> None:
        """
        Escribe un fichero escalar .scl para var_name en el timestep dado.

        Formato: %12.5e, 1 valor por línea (misma regla que la geometría;
        escribir 6/línea causaría el mismo desalineamiento de nodos).
        """
        filename = f"{self.base_name}_{var_name}_{timestep_idx:05d}.scl"
        filepath = os.path.join(self.case_dir, filename)
        with open(filepath, "w") as f:
            f.write(_truncate(f"Arritmic3D {var_name}") + "\n")
            f.write("part\n")
            f.write(f"{1:10d}\n")
            f.write("block\n")
            buf = io.StringIO()
            np.savetxt(buf, np.asarray(data, dtype=float), fmt="%12.5e")
            f.write(buf.getvalue())

    # ------------------------------------------------------------------
    # Fichero .case
    # ------------------------------------------------------------------
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
                # Variables estáticas: sin time set — un único fichero .scl
                for var_name in self._static_data:
                    filename = f"{self.base_name}_{var_name}.scl"
                    desc = _truncate(var_name, 49)
                    f.write(f"scalar per node:     {desc:<24s} {filename}\n")
                # Variables temporales: time set 1 — wildcard sobre timesteps
                for var_name in self.variables:
                    wildcard = f"{self.base_name}_{var_name}_*****.scl"
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
        """Inicializa coordenadas como índices enteros (fallback sin grid pyvista)."""
        nx = tissue.GetSizeX()
        ny = tissue.GetSizeY()
        nz = tissue.GetSizeZ()
        self._x_coords = [float(i) for i in range(nx)]
        self._y_coords = [float(i) for i in range(ny)]
        self._z_coords = [float(i) for i in range(nz)]
