from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import os
import errno
from future.moves.urllib.request import urlopen
from future.moves.urllib.error import HTTPError
from future.moves.urllib.error import URLError
from shutil import copyfileobj
from shutil import rmtree
import tarfile
import yaml
import json
import re


__all__ = ['HEPData_Manager']


class HEPData_Manager(object):

    def __init__(self, data_dir='~/.hepdata'):
        self.data_dir = _create_dir(data_dir)
        self.json_url = 'https://hepdata.net/record/{}?format=json'
        self.targz_url = 'https://www.hepdata.net/download/submission/{}/{}/{}'

    def has_data(self, inspire_id):
        return os.path.exists('{}/{}/digest.json'.format(self.data_dir,
                                                         inspire_id))

    def data_exists(self, inspire_id):
        try:
            response = urlopen(self.json_url.format(inspire_id))
            return True
        except HTTPError:
            return False
        except URLError:
            print('You are probably not connected to the internet.')
            return False

    def get_data(self, inspire_id, data_format='yaml'):
        # Checks
        #  - data exists
        #  - check for new version, if already in
        data_dir = '{}/{}'.format(self.data_dir, inspire_id)
        _create_dir(data_dir)
        digest_url = self.json_url.format(inspire_id)
        digest_file = '{}/digest.json'.format(data_dir)
        with urlopen(digest_url) as response, open(digest_file, 'wb') as out:
            copyfileobj(response, out)
        with open(digest_file) as in_file:
            data = json.load(in_file)
        targz_url = data['record']['access_urls']['links'][data_format]
        version = int(targz_url.split('/')[-2])
        targz_file = '{}/submission.tar.gz'.format(data_dir)
        with urlopen(targz_url) as response, open(targz_file, 'wb') as out:
            copyfileobj(response, out)
        tar = tarfile.open(targz_file, "r:gz")
        for x in tar.getmembers()[1:]:
            name = x.name.split('/', 1)[1]
            x_path = '{}/{}'.format(data_dir, name)
            with tar.extractfile(x) as tarf, open(x_path, 'wb') as out:
                copyfileobj(tarf, out)
        tar.close()
        for i in range(version):
            i += 1
            version_dir = '{}/v{}'.format(data_dir, i)
            _create_dir(version_dir)
            targz_url = self.targz_url.format(inspire_id, i, data_format)
            targz_file = '{}/submission.tar.gz'.format(version_dir)
            with urlopen(targz_url) as response, open(targz_file, 'wb') as out:
                copyfileobj(response, out)
            tar = tarfile.open(targz_file, "r:gz")
            for x in tar.getmembers()[1:]:
                name = x.name.split('/', 1)[1]
                x_path = '{}/{}'.format(version_dir, name)
                with tar.extractfile(x) as tarf, open(x_path, 'wb') as out:
                    copyfileobj(tarf, out)
            tar.close()

    def remove_data(self, inspire_id):
        # Checks:
        #  - data is actually stored
        data_dir = '{}/{}'.format(self.data_dir, inspire_id)
        rmtree(data_dir)

    def list_tables(self, inspire_id):
        # Checks:
        #  - data is actually stored
        digest_file = '{}/{}/digest.json'.format(self.data_dir, inspire_id)
        with open(digest_file) as in_file:
            data = json.load(in_file)
        data_tables = data['data_tables']
        return [x['processed_name'] for x in data_tables]

    def load_table(self, inspire_id, table_name, version=None,
                   data_format='yaml'):
        # Checks:
        #  - data is actually stored
        #  - Table name is valid
        #  - version is valid integer
        if version:
            table_file = '{}/{}/v{}/{}.{}'.format(self.data_dir, inspire_id,
                                                  version, table_name,
                                                  data_format)
        else:
            table_file = '{}/{}/{}.{}'.format(self.data_dir, inspire_id,
                                              table_name, data_format)
        with open(table_file) as in_file:
            data = yaml.load(in_file)
        return data


# --------------- PRIVATE --------------- #

def _create_dir(path):
    if path[0] == '~':
        path = path[2:]
        home = os.path.expanduser('~')
        path = os.path.join(home, path)
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    return str(os.path.abspath(path))
