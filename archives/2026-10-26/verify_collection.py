#!/usr/bin/env python3
"""Verify this collection without executing any archived research code."""
import hashlib
import json
import sys
import zipfile


def git_hash(kind, data):
    return hashlib.sha1(kind.encode() + b' ' + str(len(data)).encode() + b'\0' + data).hexdigest()


def tree_hash(files):
    root = {}
    for f in files:
        parts = f['path'].split('/')
        node = root
        for part in parts[:-1]:
            node = node.setdefault(part, {})
        assert parts[-1] not in node, 'Duplicate source path'
        node[parts[-1]] = (f['mode'], f['git_blob'])

    def digest(node):
        entries = []
        for name, value in node.items():
            mode, oid = ('40000', digest(value)) if isinstance(value, dict) else value
            sortkey = name.encode() + (b'/' if mode == '40000' else b'')
            entries.append((sortkey, mode.encode() + b' ' + name.encode() + b'\0' + bytes.fromhex(oid)))
        return git_hash('tree', b''.join(x[1] for x in sorted(entries)))
    return digest(root)


def verify(archive_path):
    with zipfile.ZipFile(archive_path) as z:
        names = z.namelist()
        assert len(names) == len(set(names)), 'Duplicate archive paths'
        assert z.testzip() is None, 'ZIP CRC failure'
        m = json.loads(z.read('MANIFEST.json'))
        expected = {'MANIFEST.json'}
        count = 0
        for a in m['auxiliary_files']:
            expected.add(a['path'])
            data = z.read(a['path'])
            assert len(data) == a['size_bytes']
            assert hashlib.sha256(data).hexdigest() == a['sha256'], a['path']
        for r in m['repositories']:
            commit = z.read(r['commit_object_path'])
            assert git_hash('commit', commit) == r['commit'], r['repository'] + ': commit'
            assert commit.splitlines()[0] == ('tree ' + r['root_tree']).encode()
            for f in r['files']:
                path = r['archive_prefix'] + '/' + f['path']
                expected.add(path)
                data = z.read(path)
                assert len(data) == f['size_bytes'], path + ': size'
                assert hashlib.sha256(data).hexdigest() == f['sha256'], path + ': sha256'
                assert git_hash('blob', data) == f['git_blob'], path + ': git blob'
                assert (z.getinfo(path).external_attr >> 16) == int(f['mode'], 8), path + ': mode'
                count += 1
            assert tree_hash(r['files']) == r['root_tree'], r['repository'] + ': tree'
        assert set(names) == expected, 'Unexpected or missing archive paths'
        assert count == m['summary']['tracked_files']
        with open(archive_path, 'rb') as handle:
            zip_sha256 = hashlib.file_digest(handle, 'sha256').hexdigest()
        return {'status': 'PASS', 'repositories': len(m['repositories']), 'tracked_files': count,
                'doi': m['doi'], 'zip_sha256': zip_sha256}


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: python verify_collection.py collection.zip')
    result = verify(sys.argv[1])
    print(json.dumps(result, indent=2))
