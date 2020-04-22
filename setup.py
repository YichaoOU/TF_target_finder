from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='tf_target_finder',
    version='1.0',
    author='Yichao Li',
    author_email='yl079811@ohio.edu',
    url='https://github.com/YichaoOU/TF_target_finder',
	packages=['tf_target_finder'],
    license='LICENSE',
	scripts=['assign_targets.py','co_binding_test.py','motif_scanning.py'],
	package_data={'': ["data/*"]},
	include_package_data=True,
    description=' Identify direct targets and co-binding factors',
	long_description=long_description,
	long_description_content_type='text/markdown'	,
	
)


# python setup.py sdist
# python setup.py bdist_wheel --universal
# test the distributions
# twine upload dist/*

