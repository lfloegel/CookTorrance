#ifndef PDF_H
#define PDF_H
//==============================================================================================
// Originally written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is
// distributed without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication
// along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==============================================================================================

#include "rtweekend.h"

#include "onb.h"


inline vec3 random_cosine_direction() {
    auto r1 = random_double();
    auto r2 = random_double();
    auto z = sqrt(1-r2);

    auto phi = 2*pi*r1;
    auto x = cos(phi)*sqrt(r2);
    auto y = sin(phi)*sqrt(r2);

    return vec3(x, y, z);
}


inline vec3 random_to_sphere(double radius, double distance_squared) {
    auto r1 = random_double();
    auto r2 = random_double();
    auto z = 1 + r2*(sqrt(1-radius*radius/distance_squared) - 1);

    auto phi = 2*pi*r1;
    auto x = cos(phi)*sqrt(1-z*z);
    auto y = sin(phi)*sqrt(1-z*z);

    return vec3(x, y, z);
}


class pdf  {
    public:
        virtual ~pdf() {}

        virtual double value(const vec3& direction) const = 0;
        virtual vec3 generate() const = 0;
};


class cosine_pdf : public pdf {
    public:
        cosine_pdf(const vec3& w) { uvw.build_from_w(w); }

        virtual double value(const vec3& direction) const override {
            auto cosine = dot(unit_vector(direction), uvw.w());
            return (cosine <= 0) ? 0 : cosine/pi;
        }

        virtual vec3 generate() const override {
            return uvw.local(random_cosine_direction());
        }

    public:
        onb uvw;
};


class hittable_pdf : public pdf {
    public:
        hittable_pdf(shared_ptr<hittable> p, const point3& origin) : ptr(p), o(origin) {}

        virtual double value(const vec3& direction) const override {
            return ptr->pdf_value(o, direction);
        }

        virtual vec3 generate() const override {
            return ptr->random(o);
        }

    public:
        point3 o;
        shared_ptr<hittable> ptr;
};


class mixture_pdf : public pdf {
    public:
        mixture_pdf(shared_ptr<pdf> p0, shared_ptr<pdf> p1) {
            p[0] = p0;
            p[1] = p1;
        }

        virtual double value(const vec3& direction) const override {
            return 0.5 * p[0]->value(direction) + 0.5 *p[1]->value(direction);
        }

        virtual vec3 generate() const override {
            if (random_double() < 0.5)
                return p[0]->generate();
            else
                return p[1]->generate();
        }

    public:
        shared_ptr<pdf> p[2];
};

//generate rays in all directions uniformly
class IsotropicPDF : public pdf {
    public:
	    virtual double value(const vec3& direction) const override{
            return 0.25/pi;
        }

	    virtual vec3 generate() const override{
            return random_on_unit_sphere();
        }
};

//https://computergraphics.stackexchange.com/questions/5955/where-does-the-cosine-factor-comes-from-in-the-ggx-pdf
class GGX_pdf : public pdf {
    public:
        GGX_pdf(const vec3& w, const vec3& in_dir, double r) : in(in_dir), rough(r) { 
            uvw.build_from_w(w); 
        }

        virtual double value(const vec3& direction) const override {
            vec3 h = unit_vector(direction - in);
            return GGX1(uvw.w(), h, rough) * abs(dot(h, uvw.w()));
        }

        virtual vec3 generate() const override {
            return uvw.local(random_GGX(rough));
        }

        vec3 random_GGX(double r) const {
            double r1 = random_double();
            double r2 = random_double();

            double phi = 2 * pi * r1;
            double theta = atan(r * sqrt(r2) / sqrt(1 - r2));

            double x = cos(phi) * sin(theta);
            double y = sin(phi) * sin(theta);
            double z = cos(theta);
            return vec3(x, y, z);
        }

    public:
        onb uvw;
        vec3 in;
        double rough;
};


#endif
