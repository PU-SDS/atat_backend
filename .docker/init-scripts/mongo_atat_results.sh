#!/bin/bash
mongo <<EOF
use $MONGO_INITDB_DATABASE;
db.createUser({
  user: '$MONGO_ATAT_USERNAME',
  pwd: '$MONGO_ATAT_PASSWORD',
  roles: [{
    role: 'readWrite',
    db: '$MONGO_INITDB_DATABASE'
  }]
});
db.createCollection('delete_me');
EOF