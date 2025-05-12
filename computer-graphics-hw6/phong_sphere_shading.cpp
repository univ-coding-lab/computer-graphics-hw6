//
//  sphere_scene.c
//  Rasterizer
//
//
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <cmath>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int     gNumVertices2 = 0;    // Number of 3D vertices.
int     gNumTriangles2 = 0;    // Number of triangles.
int* gIndexBuffer2 = NULL; // Vertex indices for the triangles.
unsigned char image2[512][512][3];
float z_buffer2[512][512];

struct Vec3 {
	float x, y, z;
};
Vec3* gVertexBuffer2 = nullptr;

Vec3 normalize2(const Vec3& v) {
	float len = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
	return { v.x / len, v.y / len, v.z / len };
}

float dot2(const Vec3& a, const Vec3& b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 compute_lighting2(const Vec3& pos, const Vec3& normal) {
	Vec3 ka = { 0.0f, 1.0f, 0.0f };
	Vec3 kd = { 0.0f, 0.5f, 0.0f };
	Vec3 ks = { 0.5f, 0.5f, 0.5f };
	float p = 32.0f;

	float Ia = 0.2f;
	Vec3 light_pos = { -4.0f, 4.0f, -3.0f };
	Vec3 light_color = { 1.0f, 1.0f, 1.0f };

	Vec3 L = normalize2({
		pos.x - light_pos.x,
		pos.y - light_pos.y,
		pos.z - light_pos.z
	});
	Vec3 V = normalize2({ -pos.x, -pos.y, -pos.z });
	Vec3 N = normalize2(normal);
	Vec3 R = normalize2({ L.x - 2 * dot2(N, L) * N.x, L.y - 2 * dot2(N, L) * N.y, L.z - 2 * dot2(N, L) * N.z });

	float diffuse = std::max(dot2(N, L), 0.0f);
	float specular = std::pow(std::max(dot2(R, V), 0.0f), p);

	Vec3 color = {
		Ia * ka.x + light_color.x * (kd.x * diffuse + ks.x * specular),
		Ia * ka.y + light_color.y * (kd.y * diffuse + ks.y * specular),
		Ia * ka.z + light_color.z * (kd.z * diffuse + ks.z * specular)
	};

	// Gamma correction
	color.x = std::pow(color.x, 1.0f / 2.2f);
	color.y = std::pow(color.y, 1.0f / 2.2f);
	color.z = std::pow(color.z, 1.0f / 2.2f);

	// Clamp to [0, 1]
	color.x = std::min(1.0f, std::max(0.0f, color.x));
	color.y = std::min(1.0f, std::max(0.0f, color.y));
	color.z = std::min(1.0f, std::max(0.0f, color.z));

	return color;
}

void draw_triangle_phong(
	Vec3 screen0, Vec3 screen1, Vec3 screen2,    // 스크린 좌표계 (transform_vertex 결과)
	Vec3 world0, Vec3 world1, Vec3 world2,     // 월드 좌표계 (model_transform 결과)
	Vec3 normal0, Vec3 normal1, Vec3 normal2     // 정점 normal (normalize된 vertex)
) {
	float min_x = std::min({ screen0.x, screen1.x, screen2.x });
	float max_x = std::max({ screen0.x, screen1.x, screen2.x });
	float min_y = std::min({ screen0.y, screen1.y, screen2.y });
	float max_y = std::max({ screen0.y, screen1.y, screen2.y });

	int minX = std::max(0, (int)std::floor(min_x));
	int maxX = std::min(511, (int)std::ceil(max_x));
	int minY = std::max(0, (int)std::floor(min_y));
	int maxY = std::min(511, (int)std::ceil(max_y));

	float denom = (screen1.y - screen2.y) * (screen0.x - screen2.x) +
		(screen2.x - screen1.x) * (screen0.y - screen2.y);

	for (int y = minY; y <= maxY; ++y) {
		for (int x = minX; x <= maxX; ++x) {
			float px = x + 0.5f;
			float py = y + 0.5f;

			float alpha = ((screen1.y - screen2.y) * (px - screen2.x) +
				(screen2.x - screen1.x) * (py - screen2.y)) / denom;
			float beta = ((screen2.y - screen0.y) * (px - screen2.x) +
				(screen0.x - screen2.x) * (py - screen2.y)) / denom;
			float gamma = 1.0f - alpha - beta;

			if (alpha >= 0 && beta >= 0 && gamma >= 0) {
				float z = alpha * screen0.z + beta * screen1.z + gamma * screen2.z;
				if (z < z_buffer2[y][x]) {
					z_buffer2[y][x] = z;

					// 보간된 월드 좌표
					Vec3 pos = {
						alpha * world0.x + beta * world1.x + gamma * world2.x,
						alpha * world0.y + beta * world1.y + gamma * world2.y,
						alpha * world0.z + beta * world1.z + gamma * world2.z
					};

					// 보간된 normal
					Vec3 normal = normalize2({
						alpha * normal0.x + beta * normal1.x + gamma * normal2.x,
						alpha * normal0.y + beta * normal1.y + gamma * normal2.y,
						alpha * normal0.z + beta * normal1.z + gamma * normal2.z
						});

					Vec3 color = compute_lighting2(pos, normal);
					image2[y][x][0] = (unsigned char)(color.x * 255);
					image2[y][x][1] = (unsigned char)(color.y * 255);
					image2[y][x][2] = (unsigned char)(color.z * 255);
				}
			}
		}
	}
}


Vec3 model_transform2(const Vec3& v) {
	return { v.x * 2, v.y * 2, v.z * 2 - 7 };
}

Vec3 perspective_transform2(const Vec3& v) {
	float l = -0.1f, r = 0.1f, b = -0.1f, t = 0.1f, n = -0.1f, f = -1000.0f;

	float x = (2 * n * v.x) / (r - l);
	float y = (2 * n * v.y) / (t - b);
	float z = (f + n + 2 * f * n / v.z) / (f - n);
	return { x, y, z };
}

Vec3 viewport_transform2(const Vec3& v) {
	int nx = 512, ny = 512;
	return {
		(v.x + 1.0f) * 0.5f * nx,
		(v.y + 1.0f) * 0.5f * ny,
		v.z
	};
}

Vec3 transform_vertex2(const Vec3& v) {
	Vec3 model = model_transform2(v);
	Vec3 clip = perspective_transform2(model);
	Vec3 ndc = { clip.x / -model.z, clip.y / model.z, clip.z };
	Vec3 screen = viewport_transform2(ndc);
	return screen;
}

void clear_buffers2() {
	for (int y = 0; y < 512; ++y) {
		for (int x = 0; x < 512; ++x) {
			image2[y][x][0] = image2[y][x][1] = image2[y][x][2] = 0;
			z_buffer2[y][x] = 1e9f;
		}
	}
}

void render_scene2() {
	for (int i = 0; i < gNumTriangles2; ++i) {
		int k0 = gIndexBuffer2[3 * i + 0];
		int k1 = gIndexBuffer2[3 * i + 1];
		int k2 = gIndexBuffer2[3 * i + 2];

		Vec3 world0 = model_transform2(gVertexBuffer2[k0]);
		Vec3 world1 = model_transform2(gVertexBuffer2[k1]);
		Vec3 world2 = model_transform2(gVertexBuffer2[k2]);

		Vec3 screen0 = transform_vertex2(gVertexBuffer2[k0]);
		Vec3 screen1 = transform_vertex2(gVertexBuffer2[k1]);
		Vec3 screen2 = transform_vertex2(gVertexBuffer2[k2]);

		Vec3 n0 = normalize2(gVertexBuffer2[k0]);
		Vec3 n1 = normalize2(gVertexBuffer2[k1]);
		Vec3 n2 = normalize2(gVertexBuffer2[k2]);

		draw_triangle_phong(screen0, screen1, screen2, world0, world1, world2, n0, n1, n2);
	}
}

void save_img_bmp2(const char* filename) {
	FILE* f;
	int w = 512, h = 512;
	int filesize = 54 + 3 * w * h;
	unsigned char bmpfileheader[14] = {
		'B','M', filesize & 255, (filesize >> 8) & 255, (filesize >> 16) & 255, (filesize >> 24) & 255,
		0,0, 0,0, 54,0,0,0
	};
	unsigned char bmpinfoheader[40] = {
		40,0,0,0, w & 255, (w >> 8) & 255, (w >> 16) & 255, (w >> 24) & 255,
		h & 255, (h >> 8) & 255, (h >> 16) & 255, (h >> 24) & 255,
		1,0, 24,0, 0,0,0,0,
		0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0
	};

	f = fopen(filename, "wb");
	fwrite(bmpfileheader, 1, 14, f);
	fwrite(bmpinfoheader, 1, 40, f);

	for (int y = h - 1; y >= 0; y--) {
		for (int x = 0; x < w; x++) {
			unsigned char r = image2[y][x][0];
			unsigned char g = image2[y][x][1];
			unsigned char b = image2[y][x][2];
			unsigned char color[3] = { b, g, r };
			fwrite(color, 1, 3, f);
		}
	}
	fclose(f);
}

void create_scene2()
{
	int width = 32;
	int height = 16;

	float theta, phi;
	int t;

	gNumVertices2 = (height - 2) * width + 2;
	gNumTriangles2 = (height - 2) * (width - 1) * 2;
	gVertexBuffer2 = new Vec3[gNumVertices2];

	// TODO: Allocate an array for gNumVertices2 vertices.

	gIndexBuffer2 = new int[3 * gNumTriangles2];

	t = 0;
	for (int j = 1; j < height - 1; ++j)
	{
		for (int i = 0; i < width; ++i)
		{
			theta = (float)j / (height - 1) * M_PI;
			phi = (float)i / (width - 1) * M_PI * 2;

			float   x = sinf(theta) * cosf(phi);
			float   y = cosf(theta);
			float   z = -sinf(theta) * sinf(phi);

			// TODO: Set vertex t in the vertex array to {x, y, z}.

			gVertexBuffer2[t++] = { x, y, z };
		}
	}

	// TODO: Set vertex t in the vertex array to {0, 1, 0}.

	gVertexBuffer2[t++] = { 0.0f, 1.0f, 0.0f };

	// TODO: Set vertex t in the vertex array to {0, -1, 0}.

	gVertexBuffer2[t++] = { 0.0f, -1.0f, 0.0f };

	t = 0;
	for (int j = 0; j < height - 3; ++j)
	{
		for (int i = 0; i < width - 1; ++i)
		{
			gIndexBuffer2[t++] = j * width + i;
			gIndexBuffer2[t++] = (j + 1) * width + (i + 1);
			gIndexBuffer2[t++] = j * width + (i + 1);
			gIndexBuffer2[t++] = j * width + i;
			gIndexBuffer2[t++] = (j + 1) * width + i;
			gIndexBuffer2[t++] = (j + 1) * width + (i + 1);
		}
	}
	for (int i = 0; i < width - 1; ++i)
	{
		gIndexBuffer2[t++] = (height - 2) * width;
		gIndexBuffer2[t++] = i;
		gIndexBuffer2[t++] = i + 1;
		gIndexBuffer2[t++] = (height - 2) * width + 1;
		gIndexBuffer2[t++] = (height - 3) * width + (i + 1);
		gIndexBuffer2[t++] = (height - 3) * width + i;
	}

	// The index buffer has now been generated. Here's how to use to determine
	// the vertices of a triangle. Suppose you want to determine the vertices
	// of triangle i, with 0 <= i < gNumTriangles2. Define:
	//
	// k0 = gIndexBuffer2[3*i + 0]
	// k1 = gIndexBuffer2[3*i + 1]
	// k2 = gIndexBuffer2[3*i + 2]
	//
	// Now, the vertices of triangle i are at positions k0, k1, and k2 (in that
	// order) in the vertex array (which you should allocate yourself at line
	// 27).
	//
	// Note that this assumes 0-based indexing of arrays (as used in C/C++,
	// Java, etc.) If your language uses 1-based indexing, you will have to
	// add 1 to k0, k1, and k2.
}

int main() {
	create_scene2();
	clear_buffers2();
	render_scene2();
	save_img_bmp2("phong.bmp");
	return 0;
}
