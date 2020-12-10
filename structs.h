#pragma once

struct Vertex3f { float x, y, z; }; // a single vertex position structure matching certain format

using Normal3f = Vertex3f; // a single vertex normal structure matching certain format

struct Coord2f { float u, v; }; // texture coord structure

struct Triangle3ui { unsigned int v0, v1, v2; }; // indicies of a single triangle, the struct must match certain format, e.g. RTC_FORMAT_UINT3

struct Color3f { float r, g, b; };

struct RTC_ALIGN( 16 ) Color4f {
	struct { float r, g, b, a; }; // a = 1 means that the pixel is opaque

private:
	float compress_(float u) {
		if (u <= 0) return 0;
		if (u >= 1) return 1;
		if (u <= 0.0031308f) return 12.92f * u;
		return (1.055f * pow(u, 1.0f / 2.4f) - 0.055f);
	}

	float expand_(float u) {
		if (u <= 0) return 0;
		if (u >= 1) return 1;
		if (u <= 0.04045f) return u / 12.92f;
		return pow((u + 0.055f) / 1.055f, 2.4f);
	}

	Color4f mixLinear(Color4f c0, Color4f c1) {
		float alphaDenominator = c0.a + c1.a;
		float a0 = c0.a / alphaDenominator;
		float a1 = c1.a / alphaDenominator;
		return Color4f{
			c0.r * a0 + c1.r*a1,
			c0.g * a0 + c1.g*a1,
			c0.b * a0 + c1.b*a1,
			1.0f
		};
	}

	Color4f mix_srgb(Color4f c0, Color4f c1) {
		return mixLinear(c0.expand(), c1.expand()).compress();
	}


public:
	Color4f compress() {
		return Color4f{
			compress_(r),
			compress_(g),
			compress_(b),
			a
		};
	}

	Color4f expand() {
		return Color4f{
			expand_(r),
			expand_(g),
			expand_(b),
			a
		};
	}

	Color4f operator+(const Color4f rs) {
		return mix_srgb(*this, rs);
	}

	Color4f& operator+=(const Color4f rs) {
		*this = *this + rs;
		return *this;
	}

	Color4f operator*(const Color3f rs) {
		return Color4f{
			r * rs.r,
			g * rs.g,
			b * rs.b,
			1.0f
		};
	}

};

