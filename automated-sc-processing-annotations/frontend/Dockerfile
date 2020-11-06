FROM node:10

RUN mkdir -p /usr/src/app 
WORKDIR /usr/src/app

COPY package.json /usr/src/app
COPY package-lock.json /usr/src/app 

RUN npm install
RUN npm install node-sass@4.14.1
RUN npm install react-hook-form

COPY . /usr/src/app

CMD ["npm", "start"]
