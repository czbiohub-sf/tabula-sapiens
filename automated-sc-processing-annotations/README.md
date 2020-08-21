# tabula-sapiens-as-the-reference


<img src="https://github.com/czbiohub/tabula-sapiens/blob/master/automated-sc-processing-annotations/automated-annotations.png" width="50%" height="50%">

## Steps to run:

Install Docker

```
git clone https://github.com/czbiohub/tabula-sapiens.git
cd tabula-sapiens/automated-sc-processing-annotations/
```

If you haven't downloaded the data yet:
```
./download-data.sh
```

The run: 
```
./build-and-run-image.sh
```

In your browser, navigate to:
```
localhost:8888
```
