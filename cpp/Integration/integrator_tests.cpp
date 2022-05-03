#include "integrator_tests.h"
#include "Integration.h"

void mesh_step(std::shared_ptr<Point>& p1, std::shared_ptr<Point>& p2, std::shared_ptr<Point>& p3,
                        int& Nx, int& Ny,
                        const double& dx, const double& dy,
                        int& Nxsteps, const int& Nysteps,
                        const double subdomain_dblder_limit,
                        std::map<std::pair<int, int>, const double>& evaluated_points,
                        const int& arg_ij, const double& arg_T, const int& arg_l, const int& arg_r,
                        std::function<double(int, double, double, double, int, int)> func, std::vector<Point>& points){

    Point ystep{0, dy * Nysteps};
    Point xstep{dx * Nxsteps, - dy * Nysteps};

    double f = eval_function(*p3, Nx, Ny, arg_ij, arg_T, arg_l, arg_r, func, evaluated_points);
    double f1x = eval_function(*p3 + Point{dx * abs(Nxsteps), 0}, Nx + abs(Nxsteps), Ny, arg_ij, arg_T, arg_l, arg_r, func, evaluated_points);
    double f2x = eval_function(*p3 + Point{2 * dx * abs(Nxsteps), 0}, Nx + 2 * abs(Nxsteps), Ny, arg_ij, arg_T, arg_l, arg_r, func, evaluated_points);
    double f1y = eval_function(*p3 + Point{0, dy * abs(Nysteps)}, Nx, Ny + abs(Nysteps), arg_ij, arg_T, arg_l, arg_r, func, evaluated_points);
    double f2y = eval_function(*p3 + Point{0, 2 * dy * abs(Nysteps)}, Nx, Ny + 2 * abs(Nysteps), arg_ij, arg_T, arg_l, arg_r, func, evaluated_points);
    double d2fdx2 = (f - 2 * f1x + f2x) / pow(dx * Nxsteps, 2);
    double d2fdy2 = (f - 2 * f1y + f2y) / pow(dy * Nysteps, 2);

    if ((abs(d2fdx2) + abs(d2fdy2) > subdomain_dblder_limit) && (abs(Nxsteps) > 1) && (abs(Nysteps) > 1)){ // Increase refinement
        Point sub_origin = *p3;
        int sub_Nx_origin = Nx;
        int sub_Ny_origin = Ny;
        int sub_Nxsteps = Nxsteps / 2;
        int sub_Nysteps = Nysteps / 2;
        int sub_Nx_end = Nx + Nxsteps;
        int sub_Ny_end = Ny + Nysteps;
        double sub_subdomain_dblder_limit = subdomain_dblder_limit * 2;

        mesh_adaptive(sub_origin,
                       sub_Nx_origin, sub_Ny_origin,
                       sub_Nx_end, sub_Ny_end,
                       dx, dy,
                       sub_Nxsteps, sub_Nysteps,
                       sub_subdomain_dblder_limit,
                       evaluated_points,
                       arg_ij, arg_T, arg_l, arg_r,
                       func, points);
        // Set all points to the gridpoint at the lower right corner of the subdomain that was just integrated (if Nxsteps is positive, otherwise to the lower left corner)
        //p2 = p3;
        //p1 = p2;
        //*p3 += xstep + ystep; // Only move along x-axis
        //Nx += Nxsteps;
        //eval_function(p3, Nx, Ny, arg_ij, arg_T, arg_l, arg_r, func, evaluated_points);
        //points.push_back(Point(*p3));

        p1 = std::move(p2);
        p2 = std::move(p3);
        p3 = std::shared_ptr<Point>{new Point(*p2 + ystep)};
        Ny += Nysteps;
        eval_function(p3, Nx, Ny, arg_ij, arg_T, arg_l, arg_r, func, evaluated_points);
        points.push_back(Point(*p3));

        p1 = p3;
        p2 = p3;
        *p3 += xstep; // Set all points to the point following the refined region
        Nx += Nxsteps;
        Ny -= Nysteps;
        eval_function(p3, Nx, Ny, arg_ij, arg_T, arg_l, arg_r, func, evaluated_points);
        points.push_back(Point(*p3));

    }
    else{
        p1 = std::move(p2);
        p2 = std::move(p3);
        p3 = std::shared_ptr<Point>{new Point(*p2 + ystep)};
        Ny += Nysteps;
        eval_function(p3, Nx, Ny, arg_ij, arg_T, arg_l, arg_r, func, evaluated_points);

        points.push_back(Point(*p3));

        p1 = std::move(p2);
        p2 = std::move(p3);
        p3 = std::shared_ptr<Point>{new Point(*p2 + xstep)};
        Nx += Nxsteps;
        Ny -= Nysteps;
        eval_function(p3, Nx, Ny, arg_ij, arg_T, arg_l, arg_r, func, evaluated_points);
        points.push_back(Point(*p3));
    }

}

void mesh_adaptive(const Point& origin,
                    const int& Nx_origin, const int& Ny_origin,
                    const int& Nx_end, const int& Ny_end,
                    const double& dx, const double& dy,
                    int& Nxsteps, const int& Nysteps,
                    const double& subdomain_dblder_limit,
                    std::map<std::pair<int, int>, const double>& evaluated_points,
                    const int& arg_ij, const double& arg_T, const int& arg_l, const int& arg_r,
                    std::function<double(int, double, double, double, int, int)> func, std::vector<Point>& points){

    Point ystep{0, dy * Nysteps};
    Point xstep{dx * Nxsteps, - dy * Nysteps};

    std::shared_ptr<Point> p1, p2, p3;
    p1 = std::make_shared<Point>(origin);
    p2 = p1;
    p3 = p2;
    int Nx, Ny;
    Nx = Nx_origin;
    Ny = Ny_origin;
    eval_function(p3, Nx, Ny, arg_ij, arg_T, arg_l, arg_r, func, evaluated_points);
    points.push_back(Point(*p3));

    mesh_step(p1, p2, p3, Nx, Ny, dx, dy, Nxsteps, Nysteps, subdomain_dblder_limit, evaluated_points, arg_ij, arg_T, arg_l, arg_r, func, points);
    points.push_back(Point(*p3));
    while (Ny < Ny_end){
        while (std::min(Nx_origin, Nx_end) < Nx && Nx < std::max(Nx_origin, Nx_end)){
            mesh_step(p1, p2, p3, Nx, Ny, dx, dy, Nxsteps, Nysteps, subdomain_dblder_limit, evaluated_points, arg_ij, arg_T, arg_l, arg_r, func, points);
        }
        p1 = std::move(p2);
        p2 = std::move(p3);
        p3 = std::make_shared<Point>(*p2 + ystep);
        Ny += Nysteps;
        eval_function(p3, Nx, Ny, arg_ij, arg_T, arg_l, arg_r, func, evaluated_points);
        points.push_back(Point(*p3));
        if (Ny < Ny_end){
            Nxsteps *= -1;
            mesh_step(p1, p2, p3, Nx, Ny, dx, dy, Nxsteps, Nysteps, subdomain_dblder_limit, evaluated_points, arg_ij, arg_T, arg_l, arg_r, func, points);
            points.push_back(Point(*p3));
        }
    }
}

std::vector<std::vector<double>> mesh2d(const Point& origin, const Point& end,
                                        const double& dx, const double& dy,
                                        const int& refinement_levels,
                                        const double& subdomain_dblder_limit,
                                        const int& arg_ij, const double& arg_T, const int& arg_l, const int& arg_r,
                                        std::function<double(int, double, double, double, int, int)> func){

    int Nx_origin{0}, Ny_origin{0};
    double delta_x = end.x - origin.x;
    double delta_y = end.y - origin.y;
    int Nx_end = (int) (delta_x / dx + 0.5);
    int Ny_end = (int) (delta_y / dy + 0.5);
    int Nxsteps = refinement_levels;
    int Nysteps = refinement_levels;
    std::map<std::pair<int, int>, const double> evaluated_points;
    std::vector<Point> points;
    mesh_adaptive(origin,
                 Nx_origin, Ny_origin,
                 Nx_end, Ny_end,
                 dx, dy,
                 Nxsteps, Nysteps,
                 subdomain_dblder_limit,
                 evaluated_points,
                 arg_ij, arg_T, arg_l, arg_r,
                 func, points);

    std::vector<double> x, y, z;

    for (std::vector<Point>::iterator it = points.begin(); it != points.end(); it++){
        x.push_back(it->x);
        y.push_back(it->y);
        z.push_back(it->z);
    }
    return std::vector<std::vector<double>> {x, y, z};
}

double testfun(const int ij, const double T, const double x, const double y, const int r, const int l){
    return exp(- (pow(x - 5, 2) + pow(y - 5, 2)));
}

double testfun_linear(const int ij, const double T, const double x, const double y, const int r, const int l){
    return x + y;
}

double integrator_test(double origin_x, double origin_y, double end_x, double end_y,
                       double dx, double dy, int refinement_levels, double subdomain_dblder_limit){
    int ij{1}, r{1}, l{1}; // Dummy values
    double T{1}; // Dummy values
    Point origin{origin_x, origin_y}, end{end_x, end_y};
    double val = integrate2d(origin, end, dx, dy, refinement_levels, subdomain_dblder_limit, ij, T, r, l, &testfun);
    double a = -log(testfun(ij, T, 1, 0, r, l));
    return val; // Integral of testfun on (-inf, inf) is pi / a
}

double integrator_test_linear(double origin_x, double origin_y, double end_x, double end_y,
                       double dx, double dy, int refinement_levels, double subdomain_dblder_limit){
    int ij{1}, r{1}, l{1}; // Dummy values
    double T{1}; // Dummy values
    Point origin{origin_x, origin_y}, end{end_x, end_y};
    double val = integrate2d(origin, end, dx, dy, refinement_levels, subdomain_dblder_limit, ij, T, r, l, &testfun_linear);
    return val;
}

std::vector<std::vector<double>> mesh_test(double origin_x, double origin_y, double end_x, double end_y,
                                            double dx, double dy, int refinement_levels, double subdomain_dblder_limit){
    int ij{1}, r{1}, l{1}; // Dummy values
    double T{1}; // Dummy values
    Point origin{origin_x, origin_y}, end{end_x, end_y};
    std::vector<std::vector<double>> mesh = mesh2d(origin, end, dx, dy, refinement_levels, subdomain_dblder_limit, ij, T, r, l, &testfun);
    return mesh;
}

std::vector<std::vector<double>> mesh_test_linear(double origin_x, double origin_y, double end_x, double end_y,
                       double dx, double dy, int refinement_levels, double subdomain_dblder_limit){
    int ij{1}, r{1}, l{1}; // Dummy values
    double T{1}; // Dummy values
    Point origin{origin_x, origin_y}, end{end_x, end_y};
    std::vector<std::vector<double>> mesh = mesh2d(origin, end, dx, dy, refinement_levels, subdomain_dblder_limit, ij, T, r, l, &testfun_linear);
    return mesh;
}

