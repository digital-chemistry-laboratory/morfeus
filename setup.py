import setuptools

with open("README.rst", encoding="utf8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="morfeus-ml",
    version="0.5.0",
    author="Kjell Jorner",
    author_email="kjell.jorner@gmail.com",
    description="A package to calculate molecular features.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=["morfeus"],
    python_requires=">=3.8",
    install_requires=["fire", "numpy", "scipy"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={
        "console_scripts": [
            "morfeus=morfeus.__main__:main",
        ]
    },
)
