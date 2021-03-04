#!/bin/bash
mongo <<EOF
use $RABBITMQ_MONGO_DB;
db.createUser({
  user: '$RABBITMQ_MONGO_USERNAME',
  pwd: '$RABBITMQ_MONGO_PASSWORD',
  roles: [{
    role: 'readWrite',
    db: '$RABBITMQ_MONGO_DB'
  }]
});
db.createCollection('delete_me');
EOF