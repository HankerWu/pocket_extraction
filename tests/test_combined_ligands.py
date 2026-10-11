from pathlib import Path

import pytest
from Bio.PDB import MMCIFParser, PDBParser

from pocket_extraction.extract_ligand import extract_ligand
from pocket_extraction.extract_ligand_and_pocket import extract_ligand_and_pocket


@pytest.fixture
def structure_file(tmp_path):
    lines = []
    serial = 0
    for model in (1, 2):
        lines.append(f"MODEL     {model:4d}\n")
        for record, name, residue, chain, index, x in (
            ("ATOM", "CA", "ALA", "A", 1, 0.),
            ("ATOM", "CA", "ALA", "B", 1, 20.),
            ("HETATM", "C1", "LIG", "A", 2, 1.),
            ("HETATM", "C1", "DRG", "B", 2, 21.),
        ):
            serial += 1
            lines.append(
                f"{record:<6}{serial:5d} {name:^4} {residue:3} {chain}{index:4d}    "
                f"{x:8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{20.:6.2f}           C  \n"
            )
        lines.append("ENDMDL\n")
    path = tmp_path / "input.pdb"
    path.write_text("".join(lines) + "END\n")
    return path


def read_residues(path):
    parser = PDBParser(QUIET=True) if path.suffix == ".pdb" else MMCIFParser(QUIET=True)
    structure = parser.get_structure("output", str(path))
    return [(residue.get_parent().get_parent().serial_num, residue.get_parent().id,
             residue.get_resname(), residue.id[1], len(list(residue.get_atoms())))
            for residue in structure.get_residues()]


@pytest.mark.parametrize("paired", [False, True])
@pytest.mark.parametrize("extension", ["pdb", "cif"])
@pytest.mark.parametrize("options,expected", [
    ({}, [(1, "A", "LIG", 2, 1), (1, "B", "DRG", 2, 1)]),
    ({"ignore_duplicates": False}, [(1, "A", "LIG", 2, 1), (1, "B", "DRG", 2, 1),
                                     (2, "A", "LIG", 2, 1), (2, "B", "DRG", 2, 1)]),
    ({"model_ids": [1], "chain_ids": ["B"], "ligand_names": ["DRG"]},
     [(2, "B", "DRG", 2, 1)]),
])
def test_combined_output_contains_every_selected_ligand(
    tmp_path, structure_file, paired, extension, options, expected,
):
    ligand = tmp_path / f"ligands.{extension}"
    pocket = tmp_path / f"pocket.{extension}"
    if paired:
        count = extract_ligand_and_pocket(str(structure_file), str(ligand), str(pocket),
                                         quiet=True, **options)
        pocket_residues = read_residues(pocket)
        assert all(row[2] == "ALA" for row in pocket_residues)
        assert {row[1] for row in pocket_residues} == {row[1] for row in expected}
    else:
        count = extract_ligand(str(structure_file), str(ligand), quiet=True, **options)
    assert count == len(expected)
    assert read_residues(ligand) == expected


@pytest.mark.parametrize("paired", [False, True])
def test_no_matching_ligands_does_not_publish_output(tmp_path, structure_file, paired):
    ligand, pocket = tmp_path / "ligand.pdb", tmp_path / "pocket.pdb"
    options = {"ligand_names": ["ABS"], "quiet": True}
    if paired:
        count = extract_ligand_and_pocket(str(structure_file), str(ligand), str(pocket), **options)
    else:
        count = extract_ligand(str(structure_file), str(ligand), **options)
    assert count == 0
    assert not ligand.exists()
    assert not pocket.exists()


@pytest.mark.parametrize("paired", [False, True])
@pytest.mark.parametrize("extension", ["pdb", "cif"])
def test_multi_output_retains_existing_per_ligand_files(tmp_path, structure_file, paired, extension):
    ligand, pocket = tmp_path / f"ligands.{extension}", tmp_path / f"pocket.{extension}"
    options = {"multi": True, "ignore_duplicates": False, "quiet": True}
    if paired:
        count = extract_ligand_and_pocket(str(structure_file), str(ligand), str(pocket), **options)
    else:
        count = extract_ligand(str(structure_file), str(ligand), **options)
    assert count == 4
    ligands = sorted(tmp_path.glob(f"ligands*.{extension}"))
    assert len(ligands) == 4
    assert [read_residues(path)[0][2] for path in ligands] == ["LIG", "DRG", "LIG", "DRG"]
    assert all(len(read_residues(path)) == 1 for path in ligands)
    if paired:
        pockets = sorted(tmp_path.glob(f"pocket*.{extension}"))
        assert len(pockets) == 4
        assert [{row[1] for row in read_residues(path)} for path in pockets] == [
            {"A"}, {"B"}, {"A"}, {"B"},
        ]


@pytest.mark.parametrize("paired", [False, True])
def test_actual_cli_exclude_keeps_other_ligands(tmp_path, structure_file, paired):
    import os
    import subprocess
    import sys

    module = "extract_ligand_and_pocket" if paired else "extract_ligand"
    ligand, pocket = tmp_path / "cli_ligand.pdb", tmp_path / "cli_pocket.pdb"
    arguments = [sys.executable, "-m", f"pocket_extraction.{module}", str(structure_file)]
    if paired:
        arguments += ["-ol", str(ligand), "-op", str(pocket)]
    else:
        arguments += ["-o", str(ligand)]
    arguments += ["--exclude", "LIG", "--quiet"]
    environment = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1] / "src"))
    completed = subprocess.run(arguments, env=environment, capture_output=True, text=True, timeout=10)
    assert completed.returncode == 0, completed.stderr
    assert [row[2] for row in read_residues(ligand)] == ["DRG"]
    if paired:
        assert {row[1] for row in read_residues(pocket)} == {"B"}
