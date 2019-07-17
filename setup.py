import setuptools

with open("README.md", "r") as file:
	long_description = file.read()

setuptools.setup(
	name="slime",
	version="0.0.1",
	author="Georgia Tsambos",
	author_email="gtsambos@student.unimelb.edu.au",
	description="A package that simulates tree sequences with local ancestry information.",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/tskit-dev/tskit",
	packages=setuptools.find_packages(),
	classifiers=[
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: Mac OSX and maybe others",
	],
	install_requires=["numpy", "tskit", "msprime", "pyslim"],
	setup_requires=["setuptools_scm"],
	use_scm_version={"write_to": "slime/_version.py"},
	license="MIT"
	)