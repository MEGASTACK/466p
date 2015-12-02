from PIL import Image

def main():
	images = []
	with open("jpgfilenames.txt") as filenames:
		for line in filenames:
			images.append(Image.open(line.rstrip('\n')))
	im = images[0]
	immask = im.copy()
	xmax,ymax = im.size
	for x in range(xmax):
		for y in range(ymax):
			if should_clear_pixel((x,y), im.getpixel((x,y))):
				immask.putpixel((x,y), (255,255,255))
	immask.show() # for testing

def should_clear_pixel(xy, rgb):
	x,y = xy
	r,g,b = rgb
	return (not pretty_dark(r,g,b) 
		and not greenish(r,g,b)
		and x > 15
		and y > 17)

def pretty_dark(r,g,b):
	return (30 < r + g + b < 130)

def greenish(r,g,b):
	return (65 <= r <= 105
			and 140 <= g <= 220
			and 60 <= b <= 110)

def image_path(img_filename):
	PATH = "FinalData/NP21/"
	return PATH + img_filename

if __name__ == "__main__":
	main()
