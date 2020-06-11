# tabula-hub

# Getting setup for development

## Prerequisites

1. [Node.js](https://nodejs.org/en/) v12.16.1

   _Note_: we need to use a version compatible with AWS Elastic Beanstalk for deployment

2. AWS Elastic Beanstalk CLI (to deploy to Elastic Beanstalk)

   With Python 3 -

   ```shell
   pip install awsebcli
   ```

## Setup

1. Clone the repository
2. Install Node dependencies for `backend`

   ```shell
   cd backend
   npm install
   ```

3. Install Node depdencies for `frontend`

   ```shell
   cd frontend
   npm install
   ```

## Running

To run the server, execute:

```shell
cd backend
AUTH='disabled' npm start
```

That's it - you should be able to visit http://127.0.0.1:8081

### Changing code - backend

If you modify server code, you'll have to re-run the `npm start` command above, which can get annoying. Instead, you can run the server through `nodemon`, which will watch the server files for changes and restart it for you.

```shell
npm install nodemon
nodemon --exec npm start
```

If running `nodemon` locally do

```shell
npx nodemon --exec AUTH='disabled' npm start
```

### Changing code - frontend

As you modify front-end code, it'll need to be recompiled and pushed to `backend/public/javascripts` by Webpack. By default, the `frontend` Webpack configuration will watch for changes and push updates whenever you save a file.

To start Webpack:

```shell
cd frontend
npx webpack
```
