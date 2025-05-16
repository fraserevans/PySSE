from setuptools import setup, find_packages

setup(
    name="PySSE",
    version="1.0.0",
    packages=find_packages(),
    package_data={
        "PySSE": ["pysse_module*.so"],  # Include the compiled shared library
    },
    include_package_data=True,
    install_requires=["numpy"],
    description="Python interface to the SSE (Single Star Evolution) code",
)