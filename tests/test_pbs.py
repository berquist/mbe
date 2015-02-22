import unittest

from mbe.pbs import read_nodefile
from mbe.pbs import make_nodefile_map_nodes_to_cores


class TestParseNodefile(unittest.TestCase):

    def setUp(self):
        self.nodefile_list = [
            'n179', 'n179', 'n179', 'n179',
            'n180', 'n180', 'n180', 'n180',
            'n181', 'n181', 'n181', 'n181',
            'n182', 'n182', 'n182', 'n182',
            'n183', 'n183', 'n183', 'n183',
            'n184', 'n184', 'n184', 'n184'
        ]
        self.nodefile_map = {
            'n179': 4,
            'n180': 4,
            'n181': 4,
            'n182': 4,
            'n183': 4,
            'n184': 4
        }

    def test_read_nodefile(self):
        nodefile_path = '1880021.clusman0.localdomain'
        nodefile_list = read_nodefile(nodefile_path)
        self.assertEqual(nodefile_list, self.nodefile_list)

    def test_make_nodefile_map_nodes_to_cores(self):
        nodefile_map = make_nodefile_map_nodes_to_cores(self.nodefile_list)
        self.assertEqual(nodefile_map, self.nodefile_map)


if __name__ == "__main__":
    unittest.main()
