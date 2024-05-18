from flask import Flask, render_template, request, flash

app = Flask(__name__) 
app.secret_key = "Hello"

@app.route("/")
def todo():
	return render_template("index.html")

@app.route("/PoroClicker")
def index():
	flash("what's your name?")
	return render_template("Clicker.html")

@app.route("/greet", methods=["POST", "GET"])
def greet():
	flash("Hello "+ str(request.form['name_input'] +", Great to see you"))
	return render_template("Clicker.html")
	

if __name__ == "__main__":
	app.run(debug=True)