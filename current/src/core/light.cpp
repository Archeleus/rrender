
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */


// core/light.cpp*
#include "stdafx.h"
#include "light.h"
#include "scene.h"
#include "montecarlo.h"
#include "paramset.h"
#include "sh.h"
#include "backgroundreading.h"

ImageRead im;
///test
// Light Method Definitions

Light::~Light() {
}


bool VisibilityTester::Unoccluded(const Scene *scene) const {
    return !scene->IntersectP(r);
}


Spectrum VisibilityTester::Transmittance(const Scene *scene,
        const Renderer *renderer, const Sample *sample,
        RNG &rng, MemoryArena &arena) const {
    return renderer->Transmittance(scene, RayDifferential(r), sample,
                                   rng, arena);
}


Spectrum Light::Le(const RayDifferential & ray, Sample *sample) const {
    //testing
	/*printf("here");
	Magick::Image background;
	background.read("background-picture.jpg");
	Magick::Pixels view(background);
	float i, j;
	i = 400; j = 400;*/

/*	Magick::ColorRGB t;
	t = background.pixelColor(0, 0);*/
	//TODO: Functions in ImageMagick cannot be called too many times
	//cause it uses a lot of time, especially read()
	Spectrum colour;
	float rgbv[3];
	rgbv[0] = 0.3;
	rgbv[1] = 0.4;
	rgbv[2] = 0.7;

	if((int)sample->imageX <= 800 && (int) sample->imageY <= 600){
	im.get_rgbv((int)sample->imageX, (int)sample->imageY);
	rgbv[0] = im.colours[0];
	rgbv[1] = im.colours[1];
	rgbv[2] = im.colours[2];
	}
	//printf("%d %d\n", (int)sample->imageX, (int)sample->imageY);


	return colour.FromRGB(rgbv, SPECTRUM_REFLECTANCE);

}
Spectrum Light::Le(const RayDifferential & ray) const {
	return Spectrum(0.);
		}


LightSampleOffsets::LightSampleOffsets(int count, Sample *sample) {
    nSamples = count;
    componentOffset = sample->Add1D(nSamples);
    posOffset = sample->Add2D(nSamples);
}


LightSample::LightSample(const Sample *sample,
        const LightSampleOffsets &offsets, uint32_t n) {
    Assert(n < sample->n2D[offsets.posOffset]);
    Assert(n < sample->n1D[offsets.componentOffset]);
    uPos[0] = sample->twoD[offsets.posOffset][2*n];
    uPos[1] = sample->twoD[offsets.posOffset][2*n+1];
    uComponent = sample->oneD[offsets.componentOffset][n];
    Assert(uPos[0] >= 0.f && uPos[0] < 1.f);
    Assert(uPos[1] >= 0.f && uPos[1] < 1.f);
    Assert(uComponent >= 0.f && uComponent < 1.f);
}


void Light::SHProject(const Point &p, float pEpsilon, int lmax,
        const Scene *scene, bool computeLightVisibility, float time,
        RNG &rng, Spectrum *coeffs) const {
    for (int i = 0; i < SHTerms(lmax); ++i)
        coeffs[i] = 0.f;
    uint32_t ns = RoundUpPow2(nSamples);
    uint32_t scramble1D = rng.RandomUInt();
    uint32_t scramble2D[2] = { rng.RandomUInt(), rng.RandomUInt() };
    float *Ylm = ALLOCA(float, SHTerms(lmax));
    for (uint32_t i = 0; i < ns; ++i) {
        // Compute incident radiance sample from _light_, update SH _coeffs_
        float u[2], pdf;
        Sample02(i, scramble2D, u);
        LightSample lightSample(u[0], u[1], VanDerCorput(i, scramble1D));
        Vector wi;
        VisibilityTester vis;
        Spectrum Li = Sample_L(p, pEpsilon, lightSample, time, &wi, &pdf, &vis);
        if (!Li.IsBlack() && pdf > 0.f &&
            (!computeLightVisibility || vis.Unoccluded(scene))) {
            // Add light sample contribution to MC estimate of SH coefficients
            SHEvaluate(wi, lmax, Ylm);
            for (int j = 0; j < SHTerms(lmax); ++j)
                coeffs[j] += Li * Ylm[j] / (pdf * ns);
        }
    }
}



// ShapeSet Method Definitions
ShapeSet::ShapeSet(const Reference<Shape> &s) {
    vector<Reference<Shape> > todo;
    todo.push_back(s);
    while (todo.size()) {
        Reference<Shape> sh = todo.back();
        todo.pop_back();
        if (sh->CanIntersect())
            shapes.push_back(sh);
        else
            sh->Refine(todo);
    }
    if (shapes.size() > 64)
        Warning("Area light geometry turned into %d shapes; "
            "may be very inefficient.", (int)shapes.size());

    // Compute total area of shapes in _ShapeSet_ and area CDF
    sumArea = 0.f;
    for (uint32_t i = 0; i < shapes.size(); ++i) {
        float a = shapes[i]->Area();
        areas.push_back(a);
        sumArea += a;
    }
    areaDistribution = new Distribution1D(&areas[0], areas.size());
}


ShapeSet::~ShapeSet() {
    delete areaDistribution;
}


Point ShapeSet::Sample(const Point &p, const LightSample &ls,
                       Normal *Ns) const {
    int sn = areaDistribution->SampleDiscrete(ls.uComponent, NULL);
    Point pt = shapes[sn]->Sample(p, ls.uPos[0], ls.uPos[1], Ns);
    // Find closest intersection of ray with shapes in _ShapeSet_
    Ray r(p, pt-p, 1e-3f, INFINITY);
    float rayEps, thit = 1.f;
    bool anyHit = false;
    DifferentialGeometry dg;
    for (uint32_t i = 0; i < shapes.size(); ++i)
        anyHit |= shapes[i]->Intersect(r, &thit, &rayEps, &dg);
    if (anyHit) *Ns = dg.nn;
    return r(thit);
}


Point ShapeSet::Sample(const LightSample &ls, Normal *Ns) const {
    int sn = areaDistribution->SampleDiscrete(ls.uComponent, NULL);
    return shapes[sn]->Sample(ls.uPos[0], ls.uPos[1], Ns);
}


float ShapeSet::Pdf(const Point &p, const Vector &wi) const {
    float pdf = 0.f;
    for (uint32_t i = 0; i < shapes.size(); ++i)
        pdf += areas[i] * shapes[i]->Pdf(p, wi);
    return pdf / sumArea;
}


float ShapeSet::Pdf(const Point &p) const {
    float pdf = 0.f;
    for (uint32_t i = 0; i < shapes.size(); ++i)
        pdf += areas[i] * shapes[i]->Pdf(p);
    return pdf / sumArea;
}


