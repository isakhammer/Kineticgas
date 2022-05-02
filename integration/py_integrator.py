from copy import deepcopy
import sys

if '-debug' in sys.argv or '-Debug' in sys.argv or '-d' in sys.argv:
    from integration import Integration_d as I
else:
    from integration import Integration_r as I

def getpoints(x, y, z):
    return [I.Point(x[i], y[i], z[i]) for i in range(len(x))]

def py_mesh(origin, end, dx, dy, Nxsteps, Nysteps, func, subdomainlimit, f_subdomainlimit, prev_eval_points=None):
    eval_points = {}
    if prev_eval_points is not None:
        for k, v in prev_eval_points.items():
            eval_points[k] = v

    def eval_func(x, y, Nx, Ny):
        if (Nx, Ny) in eval_points.keys():
            return eval_points[(Nx, Ny)]
        else:
            f = func(x, y)
            eval_points[(Nx, Ny)] = f
            return f

    def step(x, y, Nx, Ny):

        d2fdy2 = (eval_func(x, y + 2 * dy * Nysteps, Nx, Ny + 2 * Nysteps) - 2 * eval_func(x, y + dy * Nysteps, Nx, Ny + Nysteps) + eval_func(x, y, Nx, Ny)) / (dy * Nysteps)**2
        d2fdx2 = (eval_func(x + 2 * dx * abs(Nxsteps), y, Nx + 2 * abs(Nxsteps), Ny) - 2 * eval_func(x + dx * abs(Nxsteps), y, Nx + abs(Nxsteps), Ny) + eval_func(x, y, Nx, Ny)) / (dx * abs(Nxsteps))**2
        f = eval_func(x, y, Nx, Ny)
        if (abs(d2fdx2) + abs(d2fdy2) > subdomainlimit or f > f_subdomainlimit) and (abs(Nxsteps) > 1 and abs(Nysteps) > 1):
            sub_origin = (x, Nx, y, Ny)
            sub_Nxsteps = Nxsteps // 2
            sub_Nysteps = Nysteps // 2
            sub_end = (Nx + Nxsteps, Ny + Nysteps)
            sub_subdomainlimit = subdomainlimit * 2
            sub_f_subdomainlimit = f_subdomainlimit * 2
            subx, suby, subz, subeval = py_mesh(sub_origin,
                                          sub_end,
                                          dx, dy, sub_Nxsteps, sub_Nysteps, func,
                                          sub_subdomainlimit, sub_f_subdomainlimit,
                                          prev_eval_points=eval_points)

            for k, v in subeval.items():
                eval_points[k] = v
            xlist.extend(subx)
            ylist.extend(suby)
            zlist.extend(subz)

            x += dx * Nxsteps
            Nx += Nxsteps
            xlist.append(x)
            ylist.append(y)
            zlist.append(eval_func(x, y, Nx, Ny))
            return x, y, Nx, Ny

        y += dy * Nysteps
        Ny += Nysteps
        xlist.append(x)
        ylist.append(y)
        zlist.append(eval_func(x, y, Nx, Ny))

        y -= dy * Nysteps
        Ny -= Nysteps
        x += dx * Nxsteps
        Nx += Nxsteps
        xlist.append(x)
        ylist.append(y)
        zlist.append(eval_func(x, y, Nx, Ny))

        return x, y, Nx, Ny

    xlist = []
    ylist = []
    zlist = []
    x, origin_Nx, y, origin_Ny = origin

    Nx, Ny = origin_Nx, origin_Ny
    xend, yend = end

    x, y, Nx, Ny = step(x, y, Nx, Ny)
    while Ny < yend:

        while origin_Nx < Nx < xend or xend < Nx < origin_Nx:
            x, y, Nx, Ny = step(x, y, Nx, Ny)

        y += dy * Nysteps
        Ny += Nysteps
        xlist.append(x)
        ylist.append(y)
        zlist.append(eval_func(x, y, Nx, Ny))
        if Ny < yend:
            Nxsteps *= -1
            x, y, Nx, Ny = step(x, y, Nx, Ny)

    return deepcopy(xlist), deepcopy(ylist), deepcopy(zlist), eval_points

def py_integrate(origin, end, dx, dy, Nxsteps, Nysteps, func, subdomainlimit, f_subdomainlimit, prev_eval_points=None):
    eval_points = {}
    if prev_eval_points is not None:
        for k, v in prev_eval_points.items():
            eval_points[k] = v

    def eval_func(x, y, Nx, Ny):
        if (Nx, Ny) in eval_points.keys():
            return eval_points[(Nx, Ny)]
        else:

            f = func(x, y)
            eval_points[(Nx, Ny)] = f
            return f

    def step(x, y, Nx, Ny, integral):

        x, y = xlist[2], ylist[2]

        d2fdy2 = (eval_func(x, y + 2 * dy * Nysteps, Nx, Ny + 2 * Nysteps) - 2 * eval_func(x, y + dy * Nysteps, Nx, Ny + Nysteps) + eval_func(x, y, Nx, Ny)) / (dy * Nysteps)**2
        d2fdx2 = (eval_func(x + 2 * dx * abs(Nxsteps), y, Nx + 2 * abs(Nxsteps), Ny) - 2 * eval_func(x + dx * abs(Nxsteps), y, Nx + abs(Nxsteps), Ny) + eval_func(x, y, Nx, Ny)) / (dx * Nxsteps)**2
        f = eval_func(x, y, Nx, Ny)
        if (abs(d2fdx2) + abs(d2fdy2) > subdomainlimit or f > f_subdomainlimit) and (abs(Nxsteps) > 1 and abs(Nysteps) > 1):
            #print('Py Refining region (', x, ', ', y, '), (', Nx, ', ', Ny, '), ', abs(d2fdx2) + abs(d2fdy2), sep='')
            #print(eval_func(x, y + 2 * dy * Nysteps, Nx, Ny + 2 * Nysteps),  eval_func(x, y + dy * Nysteps, Nx, Ny + Nysteps), eval_func(x, y, Nx, Ny), (dy * Nysteps)**2)
            #print('f2x =', eval_func(x + 2 * dx * Nxsteps, y, Nx + 2 * Nxsteps, Ny),
            #      '\nf1x =', eval_func(x + dx * Nxsteps, y, Nx + Nxsteps, Ny),
            #      '\nf =', eval_func(x, y, Nx, Ny),
            #      '\ndx^2 =',  (dx * Nxsteps)**2)
            sub_origin = (x, Nx, y, Ny)
            sub_Nxsteps = Nxsteps // 2
            sub_Nysteps = Nysteps // 2
            sub_end = (Nx + Nxsteps, Ny + Nysteps)
            sub_subdomainlimit = subdomainlimit * 2
            sub_f_subdomainlimit = f_subdomainlimit * 2

            subintegral, subeval = py_integrate(sub_origin,
                                                  sub_end,
                                                  dx, dy, sub_Nxsteps, sub_Nysteps, func,
                                                  sub_subdomainlimit, sub_f_subdomainlimit,
                                                  prev_eval_points=eval_points)

            for k, v in subeval.items():
                eval_points[k] = v
            integral += subintegral

            x += dx * Nxsteps
            Nx += Nxsteps
            xlist[0] = x
            xlist[1] = x
            xlist[2] = x
            ylist[0] = y
            ylist[1] = y
            ylist[2] = y
            f = eval_func(xlist[2], ylist[2], Nx, Ny)
            zlist[0] = f
            zlist[1] = f
            zlist[2] = f
            return x, y, Nx, Ny, integral

        ylist[0] = ylist[1]
        ylist[1] = ylist[2]
        ylist[2] += dy * Nysteps
        Ny += Nysteps
        zlist[0] = zlist[1]
        zlist[1] = zlist[2]
        zlist[2] = eval_func(xlist[2], ylist[2], Nx, Ny)
        xlist[0] = xlist[1]
        xlist[1] = xlist[2]
        points = getpoints(xlist, ylist, zlist)
        integral += I.integrate_plane(points[0], points[1], points[2])

        xlist[0] = xlist[1]
        xlist[1] = xlist[2]
        xlist[2] += dx * Nxsteps
        Nx += Nxsteps

        ylist[0] = ylist[1]
        ylist[1] = ylist[2]
        ylist[2] -= dy * Nysteps
        Ny -= Nysteps

        zlist[0] = zlist[1]
        zlist[1] = zlist[2]
        zlist[2] = eval_func(xlist[2], ylist[2], Nx, Ny)

        points = getpoints(xlist, ylist, zlist)
        integral += I.integrate_plane(points[0], points[1], points[2])

        return x, y, Nx, Ny, integral

    integral = 0
    I0 = integral
    x, origin_Nx, y, origin_Ny = origin

    Nx, Ny = origin_Nx, origin_Ny
    Nx_end, Ny_end = end
    xlist = [x, x, x]
    ylist = [y, y, y]
    f = eval_func(x, y, Nx, Ny)
    zlist = [f, f, f]
    x, y, Nx, Ny, integral = step(x, y, Nx, Ny, integral)
    while Ny < Ny_end:
        I0 = integral
        while origin_Nx < Nx < Nx_end or Nx_end < Nx < origin_Nx:
            x, y, Nx, Ny, integral = step(x, y, Nx, Ny, integral)


        ylist[0] = ylist[1]
        ylist[1] = ylist[2]
        ylist[2] += dy * Nysteps
        Ny += Nysteps
        zlist[0] = zlist[1]
        zlist[1] = zlist[2]
        zlist[2] = eval_func(xlist[2], ylist[2], Nx, Ny)
        xlist[0] = xlist[1]
        xlist[1] = xlist[2]
        points = getpoints(xlist, ylist, zlist)
        integral += I.integrate_plane(points[0], points[1], points[2])
        print(xlist[2], ylist[2], integral - I0)
        if Ny < Ny_end:
            Nxsteps *= -1
            x, y, Nx, Ny, integral = step(x, y, Nx, Ny, integral)

    return integral, eval_points
