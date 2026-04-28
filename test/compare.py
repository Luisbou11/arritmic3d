"""
compare.py
==========
Verifica que la geometria generada por EnsightGoldWriter (block rectilinear)
reproduce fielmente la malla del VTK de entrada.

Uso
---
    python test/compare.py

Requiere haber ejecutado test_ventricle_ensight.py previamente para que
exista output/ventricle.case y output/ventricle.geo.
"""

import numpy as np
import pyvista as pv


VTK_REFERENCIA  = "casos/ventricle_Tagged_2.vtk"
ENSIGHT_CASE    = "output/ventricle.case"


def extract_mesh(data) -> pv.DataSet:
    """Extrae la primera malla no nula de un MultiBlock o devuelve el propio dataset."""
    if not isinstance(data, pv.MultiBlock):
        return data
    for i in range(data.n_blocks):
        block = data[i]
        if block is not None and block.n_points > 0:
            return block
    raise ValueError("El MultiBlock no contiene ninguna malla valida.")


def compare_geometries(path_vtk: str, path_case: str) -> None:
    print("=" * 60)
    print("AUDITORIA DE GEOMETRIA: VTK de entrada vs EnSight rectilinear")
    print("=" * 60)

    grid_vtk  = pv.read(path_vtk)
    grid_ens  = extract_mesh(pv.read(path_case))

    # ------------------------------------------------------------------
    # 1. Estructura
    # ------------------------------------------------------------------
    print(f"\n[ESTRUCTURA]")
    print(f"{'Propiedad':<22} {'VTK (referencia)':<22} {'EnSight (generado)':<22}")
    print("-" * 66)
    print(f"{'Tipo':<22} {type(grid_vtk).__name__:<22} {type(grid_ens).__name__:<22}")
    print(f"{'Num. puntos':<22} {grid_vtk.n_points:<22} {grid_ens.n_points:<22}")
    print(f"{'Num. celdas':<22} {grid_vtk.n_cells:<22} {grid_ens.n_cells:<22}")

    if grid_vtk.n_points != grid_ens.n_points:
        print("\nERROR CRITICO: numero de puntos distinto. Comparacion abortada.")
        return

    # ------------------------------------------------------------------
    # 2. Bounds
    # ------------------------------------------------------------------
    b_vtk = np.array(grid_vtk.bounds)
    b_ens = np.array(grid_ens.bounds)

    print(f"\n[BOUNDS]")
    labels = ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"]
    for label, v, e in zip(labels, b_vtk, b_ens):
        diff = abs(v - e)
        flag = "OK" if diff < 1e-3 else "ERROR"
        print(f"  {label:<6} VTK={v:12.5f}  ENS={e:12.5f}  diff={diff:.2e}  {flag}")

    # ------------------------------------------------------------------
    # 3. Coordenadas nodo a nodo
    # ------------------------------------------------------------------
    # PyVista devuelve los puntos de un RectilinearGrid en orden
    # x-inner, y-middle, z-outer (igual que nuestro exporter),
    # por lo que la comparacion directa es valida.
    pts_vtk = grid_vtk.points
    pts_ens = grid_ens.points

    diff     = np.abs(pts_vtk - pts_ens)
    max_diff = diff.max()
    rms_diff = np.sqrt(np.mean(diff ** 2))

    print(f"\n[COORDENADAS NODO A NODO]")
    print(f"  Diferencia maxima : {max_diff:.3e}")
    print(f"  Diferencia RMS    : {rms_diff:.3e}")

    if max_diff < 1e-3:
        print("  OK: todos los nodos coinciden dentro de la tolerancia.")
    else:
        worst = np.unravel_index(diff.max(axis=1).argmax(), (grid_vtk.n_points,))
        idx = worst[0]
        print(f"  ERROR: mayor discrepancia en nodo {idx}")
        print(f"    VTK     : {pts_vtk[idx]}")
        print(f"    EnSight : {pts_ens[idx]}")


def main() -> None:
    compare_geometries(VTK_REFERENCIA, ENSIGHT_CASE)


if __name__ == "__main__":
    main()
