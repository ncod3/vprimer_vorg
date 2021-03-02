import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="vprimer_vorg",
    version="0.0.3",
    author="satoshi-natsume",
    author_email="s-natsume@ibrc.or.jp",
    description="Vprimer_vorg",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ncod3/vprimer_vorg",
    packages=setuptools.find_packages(),
    license='GPL',
    entry_points = {
        'console_scripts': ['vprimer = vprimer.main:main']
    },
    python_requires='>=3.7',
)
