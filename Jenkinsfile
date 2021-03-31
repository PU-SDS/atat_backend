pipeline {
    agent { docker { image: 'ubuntu' } }

    stages {
        stage('Install Dependencies') {
            steps {
                echo 'Installing kubectl.'
                echo 'Updating package manager.'
                sudo apt-get update
                echo 'Installing dependencies.'
                sudo apt-get install -y apt-transport-https
                echo 'Downloading Google Cloud public signing key.'
                sudo curl -fsSLo /usr/share/keyrings/kubernetes-archive-keyring.gpg https://packages.cloud.google.com/apt/doc/apt-key.gpg
                echo 'Adding Kubernetes apt repository.'
                echo "deb [signed-by=/usr/share/keyrings/kubernetes-archive-keyring.gpg] https://apt.kubernetes.io/ kubernetes-xenial main" | sudo tee /etc/apt/sources.list.d/kubernetes.list
                echo 'Updating package index with new repository.'
                sudo apt-get update
                echo 'Installing Kubernetes.'
                sudo apt-get install -y kubectl
            }
        }
    }

}