from flask import Flask, request, render_template, jsonify, make_response, session, abort
import os
from openai import OpenAI
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import re
from io import StringIO
from flask_session import Session

# Harper Hua

app = Flask(__name__)
PASSWORD = "chatTS"
app.secret_key = os.getenv('FLASK_SECRET_KEY')

# Ensure that the secret key is set
if not app.secret_key:
    raise ValueError("No secret key set for Flask application. Set FLASK_SECRET_KEY environment variable.")


app.config['SESSION_TYPE'] = 'filesystem'
Session(app)

def create_openai_client(api_key):
    return OpenAI(api_key=api_key)

def check_api_key_validity(api_key):
    OpenAI(api_key=api_key)
    try:
        OpenAI.Model.list()
        return True
    except Exception as e:
        return False

def refine_question(client, prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            max_tokens=2000,
            temperature=0,
            messages=[
                {"role": "system", "content": "You rephrase questions to the desired format without modifying their meaning. Example: 'Does anyone have diabetes?' -> 'Does the donor have diabetes?'; 'How old are people?' -> 'What is the donor's age?'"},
                {"role": "user", "content": f"{prompt}"}
            ], 
        )
        return response.choices[0].message.content.strip()
    except Exception as e:
        print(f"An error occurred: {e}")
        return "Error: Could not generate response."

def create_chat_completion(client, prompt, file_content):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            max_tokens=2000,
            messages=[
                {"role": "system", "content": "You are a medical information retriever."},
                {"role": "user", "content": f"Here is the donor information: {file_content},"},
                {"role": "user", "content": prompt}
            ], 
        )

        return response.choices[0].message.content.strip()

    except Exception as e:
        print(f"An error occurred: {e}")
        return "Error: Could not generate response."

def create_column_name(client, prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            max_tokens=30,
            temperature=0,
            messages=[
                {"role": "system", "content": "You will be given a question. Make a concise column label. Example: 'What is the donor's age?' -> 'Age'"},
                {"role": "user", "content": f"{prompt}"}
            ], 
        )
        return response.choices[0].message.content.strip()

    except Exception as e:
        print(f"An error occurred: {e}")
        return "Error: Could not generate response."
    
def create_simplified_answer(client, answer, prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            max_tokens=200,
            messages=[
                {"role": "system", "content": f"You will be given a question and answer. Extract answer as a minimal string (only yes/no). Example prompt: 'What is the donor's age?' Answer: 'The donor is 25 years old.' Simplified answer: '25'"},
                {"role": "user", "content": f"{prompt}"},
                {"role": "user", "content": f"{answer}"}
            ], 
        )
        return response.choices[0].message.content.strip()

    except Exception as e:
        print(f"An error occurred: {e}")
        return "Error: Could not generate response."

def aggregate_donor_health_info(client, prompt, column_name, binary=False):
    donor_files_directory = '/home/ubuntu/chatTSP/row_removed/'
    donor_list = os.listdir(donor_files_directory)
    donor_numberlist = []
    donor_idlist = []
    for donor in donor_list:
        donor_id = donor.split('.')[0].split('_')[0]
        donor_idlist.append(donor_id)
        donor_number = re.findall(r'\d+', donor_id)[0]
        donor_numberlist.append(donor_number)
    donor_health_info = pd.DataFrame({'donor_id': donor_idlist}, index=donor_numberlist)
    donor_health_info[column_name] = None
    donor_health_info[f"{column_name} (short)"] = None

    def process_file(file_path):
        if file_path.endswith('.DS_Store'): 
            return None
        file_content = read_file_content(file_path)
        answer = create_chat_completion(client, prompt, file_content)
        simple_answer = create_simplified_answer(client, answer, prompt)
        return answer, simple_answer

    def process_files_in_parallel(donor_health_info, donor_files_directory):
        with ThreadPoolExecutor(max_workers=10) as executor:
            future_to_file = {executor.submit(process_file, os.path.join(donor_files_directory, i)): i for i in donor_list}

            for future in as_completed(future_to_file):
                file_name = future_to_file[future]
                try:
                    result = future.result()
                    if result is not None:
                        answer, simple_answer = result
                        donor_id = file_name.split('.')[0].split('_')[0]
                        donor_number = re.findall(r'\d+', donor_id)[0]
                        if binary:
                            donor_health_info.at[donor_number, column_name] = answer
                            donor_health_info.drop(columns=[f"{column_name} (short)"], inplace = True, errors='ignore')
                        else:
                            donor_health_info.at[donor_number, f"{column_name} (short)"] = simple_answer
                            # if "yes" in answer.lower() or "no" in answer.lower():
                            #     donor_health_info.at[donor_number, f"Simplified_answer_{column_name}"] = "Yes" if "yes" in answer.lower() else "No"
                            # else:
                            #     donor_health_info.drop(columns=[f"Binary_{column_name}"], inplace = True, errors='ignore')
                            donor_health_info.at[donor_number, column_name] = answer

                        print(f"Answer for {file_name}: {answer}")
                except Exception as exc:
                    print(f"{file_name} generated an exception: {exc}")

        return donor_health_info

    donor_health_info_df = process_files_in_parallel(donor_health_info, donor_files_directory)
    print(donor_health_info_df)
    donor_health_info.index = donor_health_info.index.astype(int)
    return donor_health_info_df

def read_file_content(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        return file.read()

@app.before_request
def check_password():
    password = request.args.get('password')
    if password != PASSWORD:
        abort(401)

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/ask', methods=['POST'])
def ask():
    question = request.form['question']
    answer_style = request.form['answerStyle']
    api_key = request.headers.get('Authorization')

    if not api_key:
        return jsonify({'error': 'API key is required.'})

    client = create_openai_client(api_key)

    #test if the API key is valid and quota is not exceeded
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            max_tokens=2000,
            messages=[
                {"role": "system", "content": "You are an AI assistant."},
                {"role": "user", "content": "This is a test question to check if the API key is valid."}
            ], 
        )
    except Exception as e:
        return jsonify({'error': str(e)})
    
    # is_valid = check_api_key_validity(api_key)
    # print("API key is valid" if is_valid else "API key is invalid")

    try:
        if answer_style == "binary":
            prompt = refine_question(client, prompt = question)
            print("MODIFIED QUESTION:", prompt)
            prompt = f"{prompt} (answer only with result in one word or number)"
            column_name = create_column_name(client, prompt = question)
            answer_pd = aggregate_donor_health_info(client, prompt, column_name, binary=True).sort_index()
        else:
            prompt = refine_question(client, prompt=question)
            print("MODIFIED QUESTION:", prompt)
            prompt = f"{prompt} if it's a binary question, answer in the format 'Yes' and give your supporting information or 'No'; if it starts with 'What' give your answers and all the matches."
            column_name = create_column_name(client, prompt)
            answer_pd = aggregate_donor_health_info(client, prompt, column_name, binary=False).sort_index()

        if 'answer_dataframe' in session:
            existing_df = pd.read_csv(StringIO(session['answer_dataframe']))
            combined_df = pd.merge(existing_df, answer_pd, on=['donor_id'], how='inner')
            unnamed_cols = [col for col in combined_df.columns if col.startswith('Unnamed')]
            combined_df = combined_df.drop(columns=unnamed_cols)
        else:
            combined_df = answer_pd
            unnamed_cols = [col for col in combined_df.columns if col.startswith('Unnamed')]
            combined_df = combined_df.drop(columns=unnamed_cols)
        
        session['answer_dataframe'] = combined_df.to_csv()
        styled_df = combined_df.replace(to_replace='\n', value='<br>', regex=True).style.set_properties(**{'text-align': 'left'}).set_table_styles(
        [{'selector': 'th', 'props': [('text-align', 'left')]},
        {'selector': 'td', 'props': [('border-top', '1px solid black')]},
        {'selector': 'th', 'props': [('border-top', '1px solid black')]},
        {'selector': 'table', 'props': [('border-collapse', 'collapse')]},
        {'selector': 'tr:last-child td', 'props': [('border-bottom', '1px solid black')]}])
        html_table = styled_df.to_html(escape=False)
        
        return jsonify({'answer': html_table})
    
    except Exception as e:
        return jsonify({'error': str(e)})

@app.route('/download')
def download():
    csv_data = session.get('answer_dataframe', 'No data available.')
    response = make_response(csv_data)
    response.headers['Content-Type'] = 'application/octet-stream'
    response.headers['Content-Disposition'] = 'attachment; filename=answer.csv'
    return response

@app.route('/reset', methods=['POST'])
def reset():
    session.pop('answer_dataframe', None)
    return jsonify({'message': 'The DataFrame has been successfully cleared.'})

@app.route('/shutdown', methods=['POST'])
def shutdown():
    print("Shutting down...")
    session.pop('answer_dataframe', None)
    session.pop('api_key', None)
    return "Shutting down..."

if __name__ == '__main__':
    app.run(debug=True)
