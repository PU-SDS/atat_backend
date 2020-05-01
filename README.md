# ATAT single

## Installation Instruction

Clone the ATAT single codes

```
$ git clone https://github.com/PU-SDS/atat-single.git
```

Build and Run Eleasticsearch and Kibana Containers

```
$ sudo sysctl -w vm.max_map_count=262144
$ cd vada-single
$ docker-compose build
$ docker-compose up -d
```

Install python supporting packages

```
$ virtualenv env
$ source env/bin/activate
$ pip install -r requirements.txt
```
