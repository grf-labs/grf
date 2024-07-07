from setuptools import setup
import os

def get_ext_modules():
    from pybind11.setup_helpers import Pybind11Extension, build_ext

    # Get a list of all .cpp files in the core/src directory and its subdirectories
    grf_src_dir = "../core/src"
    grf_sources = []
    for root, dirs, files in os.walk(grf_src_dir):
        for file in files:
            if file.endswith(".cpp"):
                grf_sources.append(os.path.join(root, file))

    ext_modules = [
        Pybind11Extension(
            "grf._grf_python",
            ["grf/_grf_python.cpp", "grf/PyUtilities.cpp"] + grf_sources,
            include_dirs=[
                "../core/src",
                "../core/third_party",
                ".",
            ],
            extra_compile_args=["-std=c++17", "-fvisibility=hidden", "-g0"],
            language="c++",
        ),
    ]
    return ext_modules

def get_cmdclass():
    from pybind11.setup_helpers import build_ext
    return {"build_ext": build_ext}

setup(
    name="grf",
    version="0.1",
    author="Apoorva Lal",
    description="Python bindings for GRF",
    ext_modules=get_ext_modules(),
    cmdclass=get_cmdclass(),
    zip_safe=False,
    packages=["grf"],
    setup_requires=[
        "pybind11",
    ]
)
