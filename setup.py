import setuptools

with open("README.md", "r") as file:
	long_description = file.read()

setuptools.setup(
	name="slime-gtsambos",
	version="0.1",
	author="Georgia Tsambos",
	author_email="gtsambos@student.unimelb.edu.au",
	description="A package that simulates tree sequences with local ancestry information.",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="tba",
	packages=setuptools.find_packages(),
	classifiers=[
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: Mac OSX and maybe others",
	]
	)