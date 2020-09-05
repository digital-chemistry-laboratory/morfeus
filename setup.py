import setuptools

with open("README.rst", "r") as fh:
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
    package_data={"morfeus": ['../data/c6_reference_data.pickle']},
    install_requires=["numpy", "scipy"],
    extras_require = {'extras': ["vtk",
                                 "pyvista",
                                 "pyvistaqt",
                                 "pymeshfix",
                                 "matplotlib"
                                 ]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={
        "console_scripts": [
            "morfeus_sterimol=morfeus.script_sterimol:main",
            "morfeus_sasa=morfeus.script_sasa:main",
            "morfeus_buried_volume=morfeus.script_buried_volume:main",
            "morfeus_cone_angle=morfeus.script_cone_angle:main",
            "morfeus_dispersion=morfeus.script_dispersion:main",
            "morfeus_local_force=morfeus.script_local_force:main",
        ]}
)
