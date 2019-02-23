import sys
from setuptools import find_packages, setup
from setuptools_rust import RustExtension


PYTHON_MAJOR_VERSION = sys.version_info[0]

setup_requires = ['setuptools-rust>=0.6.0']
install_requires = ['numpy']
test_requires = install_requires

setup(
    name='global_supertrees',
    version='0.1.0',
    description='Global supertrees extension for Python.',
    rust_extensions=[RustExtension(
        'global_supertrees.convert_file_to_distance_matrix',
        './Cargo.toml',
        rustc_flags=['--cfg=Py_{}'.format(PYTHON_MAJOR_VERSION)],
        features=['numpy/python{}'.format(PYTHON_MAJOR_VERSION)],
    )],
    install_requires=install_requires,
    setup_requires=setup_requires,
    packages=find_packages(),
    zip_safe=False,
)