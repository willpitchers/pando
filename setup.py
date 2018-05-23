'''
Uses python3.
Email: dr.mark.schultz@gmail.com
Github: schultzm
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.
You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

from setuptools import setup, find_packages
import pando
import os

LONG_DESCRIPTION = 'HAV Transmission Cluster Detection'

if os.path.exists('README'):
    LONG_DESCRIPTION = open('README').read()

setup(

    name=pando.__name__,
    version=pando.__version__,
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'pando = pando.__main__:main'
        ]
    },

    description=pando.__description__,
    long_description=LONG_DESCRIPTION,
    classifiers=['Development Status :: 3 - Alpha',
                 'License :: OSI Approved :: GNU Affero General ' +
                 'Public License v3 or later (AGPLv3+)',
                 'Programming Language :: Python :: 3.5',
                 'Topic :: Scientific/Engineering :: Bio-Informatics',
                 'Topic :: Scientific/Engineering :: Medical Science Apps.',
                 'Intended Audience :: Science/Research'],
    keywords=['Hepatitis',
              'transmission',
              'cluster'],
    download_url=pando.__install__,
    author=pando.__author__,
    author_email=pando.__author_email__,
    license=pando.__license__,
    package_data={'': ['*.fa']},
    install_requires=['ete3>=3.1.1',
                      'ruffus>=2.6.2'],
)
