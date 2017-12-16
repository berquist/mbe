"""Create a MANIFEST.in file for distributing soruces with distutils."""

from __future__ import print_function

import glob
import os.path


files = ['INSTALL.md', 'LICENSE.md', 'README.md',]
files += ['setup.py']

# source = os.path.join('src', 'cclib')
# files.append(os.path.join(source, '__init__.py'))
# files.append(os.path.join('src', 'scripts', 'ccget'))
# files.append(os.path.join('src', 'scripts', 'ccwrite'))
# files.append(os.path.join('src', 'scripts', 'cda'))

# folders = ['bridge', 'method', 'parser', 'progress', 'writer']
# for folder in folders:
#     files.extend(glob.glob(os.path.join(source, folder, '*.py')))

source = 'mbe'

for f in files:
    if not os.path.isfile(f):
        print('{} does not exist'.format(f))

with open('MANIFEST.in', 'w') as manifest_file:
    print('\n'.join(files), file=manifest_file)
