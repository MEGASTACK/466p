from PIL import Image

def main():
	images = []
	with open("jpgfilenames.txt") as filenames:
		for line in filenames:
			images.append(Image.open(line.rstrip('\n')))
	im = images[0]
	xmax,ymax = im.size
	for x in range(xmax):
		for y in range(ymax):
			if im.getpixel((x,y)) != 256:
				print im.getpixel((x,y))

def image_path(img_filename):
	PATH = "FinalData/NP21/"
	return PATH + img_filename

if __name__ == "__main__":
	main()
