import setuptools

URL = "https://github.com/digital-chemistry-laboratory/morfeus"
DESCRIPTION = "A Python package for calculating molecular features"
LONG_DESCRIPTION = f"""\
{DESCRIPTION}. For more information, see the [project repository]({URL}).
"""

setuptools.setup(
    name="morfeus-ml",
    version="0.8.0",
    author="Kjell Jorner",
    author_email="kjell.jorner@gmail.com",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url=URL,
    packages=["morfeus"],
    package_data={"morfeus": ["py.typed"]},
    python_requires=">=3.10",
    install_requires=["fire", "numpy", "packaging", "scipy"],
    license="MIT",
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
