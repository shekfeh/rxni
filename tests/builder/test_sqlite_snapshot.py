import sqlite3,tempfile,unittest
from pathlib import Path
from rxni.builder.auxiliary import _sqlite_snapshot
class SqliteSnapshotTests(unittest.TestCase):
    def test_wal_snapshot(self):
        with tempfile.TemporaryDirectory() as td:
            src=Path(td)/'source.sqlite';dst=Path(td)/'copy.sqlite'
            cx=sqlite3.connect(src)
            cx.execute('PRAGMA journal_mode=WAL')
            cx.execute('CREATE TABLE data (k INTEGER PRIMARY KEY, v TEXT)')
            cx.execute('INSERT INTO data VALUES (1,?)',('chemistry',));cx.commit()
            _sqlite_snapshot(src,dst)
            read=sqlite3.connect(dst)
            self.assertEqual(read.execute('SELECT v FROM data WHERE k=1').fetchone()[0],'chemistry')
            read.close();cx.close()
if __name__=='__main__':unittest.main()
