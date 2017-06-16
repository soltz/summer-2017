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


__all__ = ['HEPData_Manager']


class HEPData_Manager(object):

    def __init__(self, data_dir='~/.hepdata'):
        self.data_dir = _create_dir(data_dir)
        # self.data_dir = data_dir
        self.json_url = 'https://hepdata.net/record/{}?format=json'
        self.zip_url = 'https://www.hepdata.net/download/submission/{}/1/{}'

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
        digest_file = '{}/{}/digest.json'.format(self.data_dir, inspire_id)
        with urlopen(digest_url) as response, open(digest_file, 'wb') as out:
            copyfileobj(response, out)
        targz_url = self.zip_url.format(inspire_id, data_format)
        targz_file = '{}/{}/submission.tar.gz'.format(self.data_dir,
                                                      inspire_id)
        with urlopen(targz_url) as response, open(targz_file, 'wb') as out:
            copyfileobj(response, out)
        tar = tarfile.open(targz_file, "r:gz")
        for x in tar.getmembers()[1:]:
            name = x.name.split('/', 1)[1]
            x_path = '{}/{}/{}'.format(self.data_dir, inspire_id, name)
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

    def load_table(self, inspire_id, table_name, data_format='yaml'):
        # Checks:
        #  - data is actually stored
        #  - Table name is valid
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
