LDLIBS=-lpigpio -lpthread

install: noiscillate noiscillate.service
	sudo cp -v noiscillate /usr/local/sbin/
	sudo cp -v noiscillate.service /etc/systemd/system
	sudo systemctl enable noiscillate.service
	sudo systemctl restart noiscillate.service
