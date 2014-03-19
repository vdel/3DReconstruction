all:
	echo "Use 'make build' to build PhotoStereo.\nUse 'make bundler' to run Bundler.\nUse 'make stereo' to run PhotoStereo."

build:
	make -C ../../projet/src/
	
bundler:
	../../projet/bundler/RunBundler.sh
	
stereo:
	../../projet/src/PhotoStereo

clean:
	rm -f *.key.gz
	rm -f *.txt
	rm -rf bundle
	rm -rf prepare
	rm *~
