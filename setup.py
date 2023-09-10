from setuptools import setup, find_packages

setup(
    name="thztds",
    packages=find_packages(include=["thztds"]),
    version="0.1.0",
    description="Library for material analysis with THz-TDS method",
    author="Pedro Henrique Pereira Cardoso",
    author_email="pedrohpcardoso16@gmail.com",
    install_requires=["numpy==1.24.4", "pandas==2.0.3","scipy==1.10.1"],
    setup_requires=["pytest-runner==6.0.0"],
    tests_require=["pytest==7.4.0"],
    test_suite="tests",
)
