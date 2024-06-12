from setuptools import setup
from setuptools_rust import RustExtension

setup(
    name="morphing",
    version="0.1",
    rust_extensions=[RustExtension("morphing.morphing", "Cargo.toml", binding='pyo3')],
    packages=["morphing"],
    zip_safe=False,
)