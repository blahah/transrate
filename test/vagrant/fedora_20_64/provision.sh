# Setup
sudo yum install curl

# Install a multi-user RVM install with the latest ruby`
\curl -sSL https://get.rvm.io | sudo bash -s stable --ruby
source /etc/profile.d/rvm.sh

  # Add user to RVM group
sudo usermod -a -G rvm vagrant
su vagrant

# Install transrate
gem install transrate --pre

# Install dependencies
transrate --install-deps
