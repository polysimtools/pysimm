import os
import sys
import json
import inspect
import urllib
import jinja2
import webapp2

sys.path.insert(0, 'pysimm')

from pysimm import system, forcefield

JINJA_ENVIRONMENT = jinja2.Environment(
    loader=jinja2.FileSystemLoader(os.path.join(os.path.dirname(__file__), 'views')),
    extensions=['jinja2.ext.autoescape'],
    autoescape=True)

class MainPage(webapp2.RequestHandler):
    def get(self):
        template = JINJA_ENVIRONMENT.get_template('index.html')
        self.response.write(template.render({}))


class ViewerPage(webapp2.RequestHandler):
    def get(self):
        methane_json = '{"a": [{"y": 0.0, "x": 0.0, "z": 0.0, "l": "C"}, {"y": 0.7996, "x": 0.5541, "z": 0.4965, "l": "H"}, {"y": -0.8134, "x": 0.6833, "z": -0.2536, "l": "H"}, {"y": -0.3735, "x": -0.7782, "z": 0.6692, "l": "H"}, {"y": 0.3874, "x": -0.4593, "z": -0.9121, "l": "H"}], "b": [{"b": 0, "e": 1, "o": 1}, {"b": 0, "e": 2, "o": 1}, {"b": 0, "e": 3, "o": 1}, {"b": 0, "e": 4, "o": 1}]}'
        mol_json = self.request.get('json', methane_json)
        template = JINJA_ENVIRONMENT.get_template('3dviewer.html')
        self.response.write(template.render({'json': json.dumps(mol_json)}))


class SystemReadLammps(webapp2.RequestHandler):
    def post(self):
        self.response.headers.add_header("Access-Control-Allow-Origin", "*")
        input = self.request.get('input')
        s = system.read_lammps(input)
        self.response.write(s.write_yaml('string'))

class SystemWriteLammps(webapp2.RequestHandler):
    def post(self):
        self.response.headers.add_header("Access-Control-Allow-Origin", "*")
        sy = self.request.get('system')
        s = system.read_yaml(sy)
        self.response.write(s.write_lammps('string'))


class SystemWriteXyz(webapp2.RequestHandler):
    def post(self):
        self.response.headers.add_header("Access-Control-Allow-Origin", "*")
        sy = self.request.get('system')
        s = system.read_yaml(sy)
        self.response.write(s.write_xyz('string'))


class SystemReadPubchemSmiles(webapp2.RequestHandler):
    def post(self):
        self.response.headers.add_header("Access-Control-Allow-Origin", "*")
        input = self.request.get('input')
        s = system.read_pubchem_smiles(input)
        self.response.write(s.write_yaml('string'))


class SystemApplyForcefield(webapp2.RequestHandler):
    def post(self):
        self.response.headers.add_header("Access-Control-Allow-Origin", "*")
        sy = self.request.get('system')
        ff = self.request.get('forcefield')

        if not sy:
            self.response.write('cannot read system')

        if not ff:
            self.response.write('cannot read forcefield')

        s = system.read_yaml(sy)
        if ff.lower() == 'gaff':
            f = forcefield.Gaff()
        else:
            self.response.write('forcefield not supported')        

        s.apply_forcefield(f)
        self.response.write(s.write_yaml('string'))


class SystemLmps2Xyz(webapp2.RequestHandler):
    def post(self):
        self.response.headers.add_header("Access-Control-Allow-Origin", "*")
        sy = self.request.get('lammps')
        s = system.read_lammps(sy)
        self.response.write(s.write_xyz('string'))


class SystemWriteChemdoodleJson(webapp2.RequestHandler):
    def post(self):
        self.response.headers.add_header("Access-Control-Allow-Origin", "*")
        sy = self.request.get('system')
        s = system.read_yaml(sy)
        self.response.write(s.write_chemdoodle_json('string'))


class SystemPubchemSmiles2ChemdoodleJson(webapp2.RequestHandler):
    def post(self):
        self.response.headers.add_header("Access-Control-Allow-Origin", "*")
        input = self.request.get('input')
        s = system.read_pubchem_smiles(input)
        self.response.write(s.write_chemdoodle_json('string'))


app = webapp2.WSGIApplication([
    ('/', MainPage),
    ('/view', ViewerPage),
    ('/system/readlammps', SystemReadLammps),
    ('/system/readpubchemsmiles', SystemReadPubchemSmiles),
    ('/system/applyforcefield', SystemApplyForcefield),
    ('/system/writelammps', SystemWriteLammps),
    ('/system/writexyz', SystemWriteXyz),
    ('/system/lmps2xyz', SystemLmps2Xyz),
    ('/system/writechemdoodlejson', SystemWriteChemdoodleJson),
    ('/system/pubchemsmiles2chemdoodlejson', SystemPubchemSmiles2ChemdoodleJson),
], debug=True)
