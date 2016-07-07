import os
import json
import inspect
import urllib
import jinja2
import webapp2

from google.appengine.api import users

from pysimm import system, forcefield

user_list = ['test1@gmail.com', 'test2@gmail.com']

JINJA_ENVIRONMENT = jinja2.Environment(
    loader=jinja2.FileSystemLoader(os.path.dirname(__file__)),
    extensions=['jinja2.ext.autoescape'],
    autoescape=True)

class MainPage(webapp2.RequestHandler):
    def get(self):
        user = users.get_current_user()
        if user.email() in user_list:
            self.response.write('user {} is allowed here; <a href="{}">logout</a>'.format(user.email(), users.create_logout_url('/testing')))
        else:
            self.response.write('user does not have permission; <a href="{}">logout</a>'.format(users.create_logout_url('/testing')))
#        template = JINJA_ENVIRONMENT.get_template('index.html')
#        self.response.write(template.render({}))


app = webapp2.WSGIApplication([
    ('/testing', MainPage),
], debug=True)
