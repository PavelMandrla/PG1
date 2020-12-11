#include "stdafx.h"
#include "texture.h"

Texture::Texture( const char * file_name )
{
	// image format
	FREE_IMAGE_FORMAT fif = FIF_UNKNOWN;
	// pointer to the image, once loaded
	FIBITMAP * dib =  nullptr;
	// pointer to the image data
	BYTE * bits = nullptr;

	// check the file signature and deduce its format
	fif = FreeImage_GetFileType( file_name, 0 );
	// if still unknown, try to guess the file format from the file extension
	if ( fif == FIF_UNKNOWN )
	{
		fif = FreeImage_GetFIFFromFilename( file_name );
	}
	// if known
	if ( fif != FIF_UNKNOWN )
	{
		// check that the plugin has reading capabilities and load the file
		if ( FreeImage_FIFSupportsReading( fif ) )
		{
			dib = FreeImage_Load( fif, file_name );
		}
		// if the image loaded
		if ( dib )
		{
			// retrieve the image data
			bits = FreeImage_GetBits( dib );
			//FreeImage_ConvertToRawBits()
			// get the image width and height
			width_ = int( FreeImage_GetWidth( dib ) );
			height_ = int( FreeImage_GetHeight( dib ) );

			// if each of these is ok
			if ( ( bits != 0 ) && ( width_ != 0 ) && ( height_ != 0 ) )
			{				
				// texture loaded
				scan_width_ = FreeImage_GetPitch( dib ); // in bytes
				pixel_size_ = FreeImage_GetBPP( dib ) / 8; // in bytes

				data_ = new BYTE[scan_width_ * height_]; // BGR(A) format
				FreeImage_ConvertToRawBits( data_, dib, scan_width_, pixel_size_ * 8, FI_RGBA_RED_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_BLUE_MASK, TRUE );
			}

			FreeImage_Unload( dib );
			bits = nullptr;
		}
	}	
}

Texture::~Texture()
{	
	if ( data_ )
	{
		// free FreeImage's copy of the data
		delete[] data_;
		data_ = nullptr;
		
		width_ = 0;
		height_ = 0;
	}
}

Color3f Texture::valAt(int x, int y) const {
	const int offset = y * scan_width_ + x * pixel_size_;
	const float b = data_[offset] / 255.0f;
	const float g = data_[offset + 1] / 255.0f;
	const float r = data_[offset + 2] / 255.0f;
	return Color3f{ r, g, b };
}

Color3f Texture::get_texel( const float u, const float v ) const {
	
	//assert( ( u >= 0.0f && u <= 1.0f ) && ( v >= 0.0f && v <= 1.0f ) );	
	/*
	const int x = max( 0, min( width_ - 1, int( u * width_ ) ) );
	const int y = max( 0, min( height_ - 1, int( v * height_ ) ) );
	
	return valAt(x, y);
	//////////////////////
	const int offset = y * scan_width_ + x * pixel_size_;
	const float b = data_[offset] / 255.0f;
	const float g = data_[offset + 1] / 255.0f;
	const float r = data_[offset + 2] / 255.0f;
	
	return Color3f{ r, g, b };
	
	//////////////////////////////////////////////
	*/
	float x = u * float(width_);
	float y = v * float(height_);

	while (x > width_-1) x -= float(width_);
	while (x < 0) x += float(width_);
	while (y > width_-1) y -= float(width_);
	while (y < 0) y += float(height_);
	
	float x1 = floor(x);
	float x2 = ceil(x);
	x2 -= x2 >= width_ ? width_ : 0;
	float y1 = floor(y);
	float y2 = ceil(y);
	y2 = y2 >= height_ ? y1 : y2;

	auto f_xy1 = valAt(x1, y1) * float((x2 - x) / (x2 - x1)) + valAt(x2, y1) * float((x - x1) / (x2 - x1));
	auto f_xy2 = valAt(x1, y2) * ((x2 - x) / (x2 - x1)) + valAt(x2, y2) * ((x - x1) / (x2 - x1));
	return f_xy1 * ((y2 - y) / (y2 - y1)) + f_xy2 * ((y - y1) / (y2 - y1));
	

}

int Texture::width() const
{
	return width_;
}

int Texture::height() const
{
	return height_;
}
