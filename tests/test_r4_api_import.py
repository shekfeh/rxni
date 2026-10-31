"""Import-time behavior and metadata for R4; no external library is required."""
from __future__ import annotations
import importlib, json, os, subprocess, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

class TestApiImport(unittest.TestCase):
    def test_configuration_import_does_not_initialize_server(self):
        code=("import sys;sys.path.insert(0,"+repr(str(ROOT/'src'))+");"
              "import rxni.api.config as c;import rxni.api as a;"
              "assert 'fastapi' not in sys.modules;"
              "assert 'rxni.api.app' not in sys.modules;"
              "assert callable(a.create_app);"
              "assert c.ApiSettings.__name__=='ApiSettings'")
        env=os.environ.copy();env.pop('REACTION_INTELLIGENCE_LIBRARY',None);env.pop('V5',None)
        result=subprocess.run([sys.executable,'-I','-c',code],capture_output=True,text=True,env=env)
        self.assertEqual(result.returncode,0,msg=result.stderr)

    def test_api_module_import_without_library(self):
        code=("import sys;sys.path.insert(0,"+repr(str(ROOT/'src'))+");"
              "from rxni.api.app import app, create_app;"
              "assert callable(app) and callable(create_app)")
        env=os.environ.copy();env.pop('REACTION_INTELLIGENCE_LIBRARY',None);env.pop('V5',None)
        r=subprocess.run([sys.executable,'-I','-c',code],capture_output=True,text=True,env=env)
        self.assertEqual(r.returncode,0,msg=r.stderr)

    def test_pyproject_declares_extras(self):
        try:import tomllib
        except ModuleNotFoundError:import tomli as tomllib
        proj=tomllib.loads((ROOT/'pyproject.toml').read_text())['project']
        self.assertEqual(proj['version'],'0.1.0')
        self.assertTrue({'api','cli','chemkn','dev'} <= set(proj['optional-dependencies']))
        self.assertTrue({'rxni','rxni-summary','rxni-batch'} <= set(proj['scripts']))

if __name__=='__main__':unittest.main()
