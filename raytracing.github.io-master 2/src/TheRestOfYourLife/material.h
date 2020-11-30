#ifndef MATERIAL_H
#define MATERIAL_H
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

#include "pdf.h"
#include "texture.h"


struct scatter_record {
    ray specular_ray;
    bool is_specular;
    color attenuation;
    shared_ptr<pdf> pdf_ptr;
};


class material {
    public:
        virtual color emitted(
            const ray& r_in, const hit_record& rec, double u, double v, const point3& p
        ) const {
            return color(0,0,0);
        }

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, scatter_record& srec
        ) const {
            return false;
        }

        virtual double scattering_pdf(
            const ray& r_in, const hit_record& rec, const ray& scattered
        ) const {
            return 0;
        }
};


class lambertian : public material {
    public:
        lambertian(const color& a) : albedo(make_shared<solid_color>(a)) {}
        lambertian(shared_ptr<texture> a) : albedo(a) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, scatter_record& srec
        ) const override {
            srec.is_specular = false;
            srec.attenuation = albedo->value(rec.u, rec.v, rec.p);
            srec.pdf_ptr = make_shared<cosine_pdf>(rec.normal);
            return true;
        }

        double scattering_pdf(
            const ray& r_in, const hit_record& rec, const ray& scattered
        ) const override {
            auto cosine = dot(rec.normal, unit_vector(scattered.direction()));
            return cosine < 0 ? 0 : cosine/pi;
        }

    public:
        shared_ptr<texture> albedo;
};


class metal : public material {
    public:
        metal(const color& a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, scatter_record& srec
        ) const override {
            vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
            srec.specular_ray =
                ray(rec.p, reflected + fuzz*random_in_unit_sphere(), r_in.time());
            srec.attenuation = albedo;
            srec.is_specular = true;
            srec.pdf_ptr = nullptr;
            return true;
        }

    public:
        color albedo;
        double fuzz;
};


class dielectric : public material {
    public:
        dielectric(double index_of_refraction) : ir(index_of_refraction) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, scatter_record& srec
        ) const override {
            srec.is_specular = true;
            srec.pdf_ptr = nullptr;
            srec.attenuation = color(1.0, 1.0, 1.0);
            double refraction_ratio = rec.front_face ? (1.0/ir) : ir;

            vec3 unit_direction = unit_vector(r_in.direction());
            double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
            double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

            bool cannot_refract = refraction_ratio * sin_theta > 1.0;
            vec3 direction;

            if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
                direction = reflect(unit_direction, rec.normal);
            else
                direction = refract(unit_direction, rec.normal, refraction_ratio);

            srec.specular_ray = ray(rec.p, direction, r_in.time());
            return true;
        }

    public:
        double ir; // Index of Refraction

    private:
        static double reflectance(double cosine, double ref_idx) {
            // Use Schlick's approximation for reflectance.
            auto r0 = (1-ref_idx) / (1+ref_idx);
            r0 = r0*r0;
            return r0 + (1-r0)*pow((1 - cosine),5);
        }
};


class diffuse_light : public material {
    public:
        diffuse_light(shared_ptr<texture> a) : emit(a) {}
        diffuse_light(color c) : emit(make_shared<solid_color>(c)) {}

        virtual color emitted(
            const ray& r_in, const hit_record& rec, double u, double v, const point3& p
        ) const override {
            if (!rec.front_face)
                return color(0,0,0);
            return emit->value(u, v, p);
        }

    public:
        shared_ptr<texture> emit;
};


class isotropic : public material {
    public:
        isotropic(color c) : albedo(make_shared<solid_color>(c)) {}
        isotropic(shared_ptr<texture> a) : albedo(a) {}

        //#if 0
        // Issue #669
        // This method doesn't match the signature in the base `material` class, so this one's
        // never actually called. Disabling this definition until we sort this out.

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, scatter_record& srec
        ) const override {
            srec.pdf_ptr = make_shared<IsotropicPDF>();
            srec.is_specular = false;
            srec.specular_ray = ray(rec.p, random_on_unit_sphere());
            srec.attenuation = albedo->value(rec.u, rec.v, rec.p);
            return true;
        }

        double scattering_pdf(
            const ray& r_in, const hit_record& rec, const ray& scattered
        ) const override {
            return 0.25/pi;
        }
       // #endif

    public:
        shared_ptr<texture> albedo;
};

//http://www.codinglabs.net/article_physically_based_rendering_cook_torrance.aspx
//https://computergraphics.stackexchange.com/questions/4394/path-tracing-the-cook-torrance-brdf?rq=1
class CookTorrance : public material {

    public:
        //constructor
        CookTorrance(double r, double m, double ior, const color& a)
            : rough((r = r + 0.001) > 1 ? 1 : r), metallic(m), ir(ior), albedo(a) {};

        //scattered pdf
        double scattering_pdf(
            const ray& r_in, const hit_record& rec, const ray& scattered
        ) const override {
            auto cosine = dot(rec.normal, unit_vector(scattered.direction()));
            return cosine < 0 ? 0 : (cosine/pi); //if in hemisphere
        }
            
        //scatter 
        virtual bool scatter(
            const ray& r_in, const hit_record& rec, scatter_record& srec
        ) const override {
            //outgoing direction
            vec3 wo = unit_vector(r_in.direction());

            //GGX pdf
            GGX_pdf p(rec.normal, r_in.direction(), rough);
            srec.specular_ray = ray(rec.p, p.generate());

            //incoming direction
            vec3 wi = unit_vector(srec.specular_ray.direction());
            vec3 half_vector = unit_vector(wi + wo);
            double cosine = dot(wi, rec.normal);

            //fresnel schlick
            vec3 f0 = color_normal(ir, albedo, metallic);
            f0 = fresnel_schlick(cosine, f0);

            vec3 unit(1, 1, 1);

            //F
            auto F = f0;

            //D
            auto D = GGX1(rec.normal, half_vector, rough);

            //G
            auto G = GGX_partialgeomtry(rec.normal, half_vector, wi, wo);

            //specular
            auto specular = F * (D * G / (4 * abs(dot(rec.normal, wo)) 
            * abs(dot(rec.normal, wi)) * p.value(srec.specular_ray.direction())));

            //diffuse
            auto diffuse = unit - f0 * albedo * (1/pi);

            srec.attenuation = specular + diffuse * clamp(dot(rec.normal, wi), 0);
            srec.is_specular = 1;
            srec.pdf_ptr = nullptr;
            return true;
        }
        
    public:
        //parameters
        double rough, metallic, ir;
        color albedo;

    private:
        //geometry cook torrance
        double geometry(const vec3& n, const vec3& h, const vec3& wi, const vec3& wo) const {
            double n_wi = dot(n, wi);
            double n_wo = dot(n, wo);
            double x = 2 * dot(n, h) / dot(wo, h);
            return fmin(1, fmin(x * n_wo, x * n_wi));
        }

        //GGX geometry: attenuation of light
        double GGX_partialgeomtry(const vec3& n, const vec3& h, const vec3& wi, const vec3& wo) const {
            double wo2 = dot(wo, h);
            double chi = clamp(wo2 / dot(wo, n), 0, 1);
            wo2 = wo2 * wo2;
            return (chi * 2) / (1 + sqrt(1 + rough * rough * ((1 - wo2) / wo2)));
        }


        static vec3 color_normal(double ref_idx, const color& a, double m) {
            double e = (1 - ref_idx) / (1 + ref_idx);
            vec3 f0(e, e, e);
            f0 = f0 * f0;
            return lerp(f0, a, m);;
        }

        static vec3 fresnel_schlick(double cosine, const vec3& f0) {
            vec3 x(1,1,1);
            return f0 + (x - f0) * pow(1 - cosine, 5);
        }
};


#endif
