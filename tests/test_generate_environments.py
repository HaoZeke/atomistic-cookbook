"""Unit tests for src/generate_environments.py."""

from __future__ import annotations

import io
import sys
import tempfile
import unittest
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.generate_environments import (  # noqa: E402
    CONDA_CPU_EXTRAS,
    CPU_INDEX_FLAG,
    PYTORCH_CPU_INDEX,
    append_environment_downloads,
    conda_package_name,
    main,
    render_cpu_environment,
    write_environment_variants,
)


EON_ENV = ROOT / "examples" / "eon-pet-neb" / "environment.yml"

EON_LIKE = """\
channels:
  - metatensor
  - conda-forge
  - https://repo.prefix.dev/rg-forge
dependencies:
  - ira
  - python=3.13
  - pip
  - uv
  - pip:
    - "pyeonclient[models]>=0.3.5,<0.4"
    - torch>=2.10,<2.14
    - atomistic-cookbook-utils >=0.1,<0.2
"""

ALREADY_PINNED = """\
channels:
  - conda-forge
dependencies:
  - python=3.12
  - libtorch
  - pytorch-cpu
  - lammps-metatomic=2026.07.04.mta1
  - lammps-metatomic * cpu*nompi*
  - pip
  - pip:
    - --extra-index-url https://download.pytorch.org/whl/cpu
    - torch>=2.12,<2.13
variables:
  PIP_EXTRA_INDEX_URL: https://download.pytorch.org/whl/cpu
"""

CONDA_ENGINES = """\
channels:
  - metatensor
  - conda-forge
dependencies:
  - python=3.12
  - lammps-metatomic==2026.03.30.mta1
  - plumed-metatomic=2.10.1.dev1
  - pip
  - pip:
    - ase
"""


class CondaNameTests(unittest.TestCase):
    def test_plain_and_versioned(self):
        self.assertEqual(conda_package_name("libtorch"), "libtorch")
        self.assertEqual(conda_package_name("python=3.13"), "python")
        self.assertEqual(
            conda_package_name("lammps-metatomic==2026.03.30.mta1"),
            "lammps-metatomic",
        )

    def test_channel_and_build_string(self):
        self.assertEqual(
            conda_package_name("metatensor::lammps-metatomic=*=cpu*"),
            "lammps-metatomic",
        )
        self.assertEqual(
            conda_package_name("lammps-metatomic * cpu*nompi*"),
            "lammps-metatomic",
        )


class RenderCpuTests(unittest.TestCase):
    def test_adds_pip_index_and_variable(self):
        env = yaml.safe_load(EON_LIKE)
        cpu = render_cpu_environment(env)
        pip_items = next(
            dep["pip"] for dep in cpu["dependencies"] if isinstance(dep, dict)
        )
        self.assertEqual(pip_items[0], CPU_INDEX_FLAG)
        self.assertIn("torch>=2.10,<2.14", pip_items)
        self.assertEqual(cpu["variables"]["PIP_EXTRA_INDEX_URL"], PYTORCH_CPU_INDEX)

    def test_does_not_edit_the_input_mapping(self):
        env = yaml.safe_load(EON_LIKE)
        render_cpu_environment(env)
        pip_items = next(
            dep["pip"] for dep in env["dependencies"] if isinstance(dep, dict)
        )
        self.assertNotIn(CPU_INDEX_FLAG, pip_items)
        self.assertNotIn("variables", env)

    def test_idempotent_when_already_pinned(self):
        env = yaml.safe_load(ALREADY_PINNED)
        cpu = render_cpu_environment(env)
        pip_items = next(
            dep["pip"] for dep in cpu["dependencies"] if isinstance(dep, dict)
        )
        self.assertEqual(pip_items.count(CPU_INDEX_FLAG), 1)
        conda_deps = [dep for dep in cpu["dependencies"] if isinstance(dep, str)]
        self.assertEqual(conda_deps.count("pytorch-cpu"), 1)
        self.assertEqual(conda_deps.count(CONDA_CPU_EXTRAS["lammps-metatomic"]), 1)

    def test_adds_known_conda_cpu_extras(self):
        env = yaml.safe_load(CONDA_ENGINES)
        cpu = render_cpu_environment(env)
        conda_deps = [dep for dep in cpu["dependencies"] if isinstance(dep, str)]
        self.assertIn(CONDA_CPU_EXTRAS["lammps-metatomic"], conda_deps)
        self.assertIn(CONDA_CPU_EXTRAS["plumed-metatomic"], conda_deps)
        pip_index = next(
            i
            for i, dep in enumerate(cpu["dependencies"])
            if isinstance(dep, dict) and "pip" in dep
        )
        extra_indexes = [
            i
            for i, dep in enumerate(cpu["dependencies"])
            if dep in CONDA_CPU_EXTRAS.values()
        ]
        self.assertTrue(all(index < pip_index for index in extra_indexes))

    def test_no_pip_section_still_adds_conda_extras(self):
        env = {
            "channels": ["conda-forge"],
            "dependencies": ["python=3.12", "libtorch"],
        }
        cpu = render_cpu_environment(env)
        self.assertIn("pytorch-cpu", cpu["dependencies"])
        self.assertNotIn("variables", cpu)

    def test_empty_mapping(self):
        self.assertEqual(render_cpu_environment({}), {})


class WriteVariantsTests(unittest.TestCase):
    def test_default_is_byte_copy_cpu_is_rendered(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        source = Path(tmp.name) / "environment.yml"
        source.write_text(EON_LIKE, encoding="utf-8")
        dest_dir = Path(tmp.name) / "out"

        paths = write_environment_variants(source, dest_dir)

        self.assertEqual(paths["default"].read_bytes(), source.read_bytes())
        cpu_text = paths["cpu"].read_text(encoding="utf-8")
        self.assertIn("generate_environments.py", cpu_text)
        loaded = yaml.safe_load(cpu_text)
        pip_items = next(
            dep["pip"] for dep in loaded["dependencies"] if isinstance(dep, dict)
        )
        self.assertEqual(pip_items[0], CPU_INDEX_FLAG)

    def test_eon_pet_neb_is_a_consumer(self):
        self.assertTrue(EON_ENV.is_file(), EON_ENV)
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        paths = write_environment_variants(EON_ENV, tmp.name)
        source = yaml.safe_load(EON_ENV.read_text(encoding="utf-8"))
        cpu = yaml.safe_load(paths["cpu"].read_text(encoding="utf-8"))
        source_pip = next(
            dep["pip"] for dep in source["dependencies"] if isinstance(dep, dict)
        )
        cpu_pip = next(
            dep["pip"] for dep in cpu["dependencies"] if isinstance(dep, dict)
        )
        self.assertNotIn(CPU_INDEX_FLAG, source_pip)
        self.assertEqual(cpu_pip[0], CPU_INDEX_FLAG)
        self.assertIn("torch>=2.10,<2.14", cpu_pip)


class DownloadRstTests(unittest.TestCase):
    def test_appends_once(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        rst = Path(tmp.name) / "recipe.rst"
        rst.write_text("Title\n=====\n\nBody.\n", encoding="utf-8")
        append_environment_downloads(rst)
        append_environment_downloads(rst)
        text = rst.read_text(encoding="utf-8")
        self.assertEqual(text.count("environment-cpu.yml"), 2)
        self.assertIn("environment-default.yml", text)


class CliTests(unittest.TestCase):
    def test_main_writes_both_files(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        source = Path(tmp.name) / "environment.yml"
        source.write_text(EON_LIKE, encoding="utf-8")
        out = Path(tmp.name) / "out"
        buf = io.StringIO()
        old = sys.stdout
        sys.stdout = buf
        try:
            status = main([str(source), "-o", str(out)])
        finally:
            sys.stdout = old
        self.assertEqual(status, 0)
        self.assertTrue((out / "environment-cpu.yml").is_file())
        self.assertTrue((out / "environment-default.yml").is_file())
        self.assertIn("cpu:", buf.getvalue())


if __name__ == "__main__":
    unittest.main()
