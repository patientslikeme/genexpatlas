from setuptools import setup, find_packages

setup(name='genexpatlas',
      version='0.6.0-alpha',
      description='Python interface to the Genetic Expression Atlas',
      author='Kyle Stratis',
      author_email='kstratis@patientslikeme.com',
      packages=find_packages(),
      install_requires=[
          'requests',
          'pandas',
          'xmltodict',
      ],
      project_urls={
          'PatientsLikeMe': 'https://blog.patientslikeme.com',
          'Bug reports': 'https://github.com/patientslikeme/genexpatlas/issues',
          'Source': 'https://github.com/patientslikeme/genexpatlas',
      },
)
