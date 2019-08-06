from setuptools import setup,find_packages

config = {
    'include_package_data': True,
    'description': 'SAM file alignment statistics at the read level',
    'version': '0.2.1',
    'packages': ['SAMstats','SAMstatsParallel'],
    'setup_requires': [],
    'install_requires': ['multiprocess'],
    'dependency_links': [],
    'scripts': [],
    'entry_points': {'console_scripts': ['SAMstats = SAMstats.__init__:main',
                                         'SAMstatsParallel= SAMstatsParallel.__init__:main']},
    'name': 'SAMstats'
}

if __name__== '__main__':
    setup(**config)
