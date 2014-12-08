su vagrant

# Setup
sudo apt-get update -y
sudo apt-get install -y curl build-essential git zlib1g zlib1g-dev cmake

# Install a multi-user RVM install with the latest ruby
sudo gpg --keyserver hkp://keys.gnupg.net --recv-keys D39DC0E3
\curl -sSL https://get.rvm.io | bash -s stable --ruby
source /usr/local/rvm/scripts/rvm

# Add user to RVM group
sudo usermod -a -G rvm vagrant

# Install transrate
git clone https://github.com/Blahah/transrate.git
cd transrate
sudo chmod -R 777 .

# Install transrate-tools
cd ..
git clone https://github.com/Blahah/transrate-tools.git

# Install dependencies
gem install bundler rake
bundle install
bundle exec rake compile
bundle exec bin/transrate --install-deps
