from flask import Flask
import time

app = Flask(__name__)
@app.route('/time')
def get_current_time():
    return {'time': time.time()}

@app.route('/')
def hello_world():
    return "SUH"

#if __name__ = "__main__":
#    app.run()
