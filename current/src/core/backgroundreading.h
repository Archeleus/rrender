/*
 * backgroundreading.h
 *
 *  Created on: 21 Sep 2010
 *      Author: alex
 */

#ifndef BACKGROUNDREADING_H_
#define BACKGROUNDREADING_H_

#include <ImageMagick/Magick++.h>

class ImageRead{
public:
	ImageRead();
	void get_rgbv(int, int);
public:
	float colours[3];
	Magick::Image backgroundd;
	//Magick::Pixels view(Magick::Image);
	Magick::ColorRGB t;
};

ImageRead::ImageRead(){
	backgroundd.magick("jpg");
	backgroundd.read("background-picture.jpg");
	//view(backgroundd);

}

void ImageRead::get_rgbv(int x, int y){
	Magick::Pixels view(backgroundd);

	t = backgroundd.pixelColor(x, y);
	colours[0] = t.red();
	colours[1] = t.green();
	colours[2] = t.blue();

}
#endif /* BACKGROUNDREADING_H_ */
