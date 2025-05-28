// Compile with: g++ -O2 -o a.out  main.cpp -lfftw3 -lfftw3_threads -lm

#include "./customize/params.cpp"
#include "./utils/utils.cpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <vector>
#include <string>
#include <bits/stdc++.h>


// Function to initialize wave NUmbers with FFTW ordering.
void init_wave_numbers(std::vector<double>& kx, std::vector<double>& ky) {
    kx.resize(GRID_X);
    ky.resize(GRID_Y/2+1);  
    const double dkx = 2.0 * M_PI / LENGTH_X;
    const double dky = 2.0 * M_PI / LENGTH_Y;
    for (int i = 0; i < GRID_X; ++i) {
        if (i <= GRID_X/2)
            kx[i] = i * dkx;
        else
            kx[i] = (i - GRID_X) * dkx;
    }
    for (int j = 0; j < GRID_Y/2+1; ++j) {
        ky[j] = j * dky;
    }
}

int main()
{
    //Init Threads
    fftw_init_threads();
    fftw_plan_with_nthreads(1);

    // Allocate arrays in physical (real) space.
    int real_size = GRID_X * GRID_Y;
    double* omega   = fftw_alloc_real(real_size); 
    double* N_phys  = fftw_alloc_real(real_size);
    double* A_phys  = fftw_alloc_real(real_size);

    // Allocate arrays for velocity components (physical space).
    double* u_x_phys = fftw_alloc_real(real_size);
    double* u_y_phys = fftw_alloc_real(real_size);

    // Allocate arrays in Fourier (complex) space.
    // The complex arrays are of size: GRID_X x (GRID_Y/2+1)
    int complex_size = GRID_X * (GRID_Y/2 + 1);
    fftw_complex* omega_hat = fftw_alloc_complex(complex_size);
    fftw_complex* A_hat     = fftw_alloc_complex(complex_size);
    fftw_complex* N_hat     = fftw_alloc_complex(complex_size);

    // Allocate arrays for FFTs of velocity components (for computing initial vorticity).
    fftw_complex* ux_hat = fftw_alloc_complex(complex_size);
    fftw_complex* uy_hat = fftw_alloc_complex(complex_size);

    // Plans for forward and inverse FFTs.
    fftw_plan plan_forward_omega   = fftw_plan_dft_r2c_2d(GRID_X, GRID_Y, omega, omega_hat, FFTW_MEASURE);
    fftw_plan plan_backward_omega  = fftw_plan_dft_c2r_2d(GRID_X, GRID_Y, omega_hat, omega, FFTW_MEASURE);
    fftw_plan plan_backward_A = fftw_plan_dft_c2r_2d(GRID_X, GRID_Y, A_hat, A_phys, FFTW_MEASURE);
    fftw_plan plan_backward_ux = fftw_plan_dft_c2r_2d(GRID_X, GRID_Y, ux_hat, u_x_phys, FFTW_MEASURE);
    fftw_plan plan_backward_uy = fftw_plan_dft_c2r_2d(GRID_X, GRID_Y, uy_hat, u_y_phys, FFTW_MEASURE);

    // Prepare wave NUmber arrays.
    std::vector<double> kx, ky;
    init_wave_numbers(kx, ky);
    double ens = 0;
    double ens_old = 0;

    // Output file
    std::ofstream file_velocity;
    std::ofstream file_velocity_max;
    std::ofstream file_ens;
    file_velocity.open("./output/velocity.txt");
    file_velocity_max.open("./output/velocity_max.txt");
    file_ens.open("./output/enstrophie.txt");
    
    // ------------------------------
    // INITIALIZATION FROM VELOCITY FIELD
    // ------------------------------
    for (int i = 0; i < GRID_X; ++i) {
        double x = i * LENGTH_X / GRID_X;
        for (int j = 0; j < GRID_Y; ++j) {
            double y = j * LENGTH_Y / GRID_Y;
            int idx = i*GRID_Y + j;
            double rand_vel = ((rand() % 101)-50)*AMP/50.0;
            if (j < (GRID_Y / 2))
            {
                u_x_phys[idx] = 1.0;
                u_y_phys[idx] = rand_vel;
            }
            else {
                u_x_phys[idx] = -1.0;
                u_y_phys[idx] = rand_vel;
            }
        }
    }
    // Create forward FFT plans for the velocity components.
    fftw_plan plan_forward_ux = fftw_plan_dft_r2c_2d(GRID_X, GRID_Y, u_x_phys, ux_hat, FFTW_MEASURE);
    fftw_plan plan_forward_uy = fftw_plan_dft_r2c_2d(GRID_X, GRID_Y, u_y_phys, uy_hat, FFTW_MEASURE);
    fftw_execute(plan_forward_ux);
    fftw_execute(plan_forward_uy);
    fftw_destroy_plan(plan_forward_ux);
    fftw_destroy_plan(plan_forward_uy);

    // Compute initial vorticity in Fourier space via:
    //      ω̂(k) = i*k_x * (u_y)_hat(k) - i*k_y * (u_x)_hat(k)
    // Note: Multiplication by i transforms a complex NUmber f = a + i b into i f = -b + i a.
    for (int i = 0; i < GRID_X; ++i) {
        for (int j = 0; j < GRID_Y/2+1; ++j) {
            int idx = i * (GRID_Y/2+1) + j;
            double kx_val = kx[i];
            double ky_val = ky[j];
            double a = ux_hat[idx][0], b = ux_hat[idx][1]; // Fourier coeff of u_x
            double c = uy_hat[idx][0], d = uy_hat[idx][1]; // Fourier coeff of u_y
            // Compute:
            //   i*kx*(u_y)_hat = -kx * d + i*kx * c
            //   i*ky*(u_x)_hat = -ky * b + i*ky * a
            // Therefore, omega_hat = i*kx*(u_y)_hat - i*ky*(u_x)_hat:
            omega_hat[idx][0] = (-kx_val * d) - (-ky_val * b);  // real part = -kx*d + ky*b
            omega_hat[idx][1] = (kx_val * c) - (ky_val * a);      // imag part = kx*c - ky*a
        }
    }
    // Inverse FFT to obtain the physical-space initial vorticity.
    fftw_execute(plan_backward_omega);
    for (int i = 0; i < real_size; ++i)
        omega[i] /= (GRID_X * GRID_Y);
    // At this point, omega (physical) and omega_hat (spectral) contain the initial condition.
    // ------------------------------
    // END INITIALIZATION
    // ------------------------------

    // Main time-stepping loop.
    for (int step = 0; step < STEPS_TOTAL; ++step) {
        ens_old = ens;

        // --- Step 1: Forward FFT of omega (update omega_hat from omega) ---
        fftw_execute(plan_forward_omega);

        // --- Step 2: Compute A_hat from omega_hat: A_hat = -omega_hat/|k|^2 (with care at k=0) ---
        for (int i = 0; i < GRID_X; ++i) {
            for (int j = 0; j < GRID_Y/2+1; ++j) {
                int idx = i*(GRID_Y/2+1) + j;
                double k2 = kx[i]*kx[i] + ky[j]*ky[j];
                if (k2 < 1e-10) {
                    A_hat[idx][0] = 0.0;
                    A_hat[idx][1] = 0.0;
                } else {
                    A_hat[idx][0] = -omega_hat[idx][0] / k2;
                    A_hat[idx][1] = -omega_hat[idx][1] / k2;
                }
            }
        }

        // --- Step 3: Compute velocity components from A_hat in Fourier space ---
        // u_x = ∂A/∂y  ->  Fourier: û_x = i*ky*A_hat,
        // u_y = -∂A/∂x ->  Fourier: û_y = -i*kx*A_hat.
        for (int i = 0; i < GRID_X; ++i) {
            for (int j = 0; j < GRID_Y/2+1; ++j) {
                int idx = i*(GRID_Y/2+1) + j;
                double reA = A_hat[idx][0];
                double imA = A_hat[idx][1];
                // u_x Fourier component:
                ux_hat[idx][0] = -ky[j]*imA;
                ux_hat[idx][1] =  ky[j]*reA;
                // u_y Fourier component:
                uy_hat[idx][0] =  kx[i]*imA;
                uy_hat[idx][1] = -kx[i]*reA;
            }
        }
        // Inverse FFTs to obtain u_x and u_y in physical space.
        fftw_execute(plan_backward_ux);
        fftw_execute(plan_backward_uy);
        for (int i = 0; i < real_size; ++i) {
            u_x_phys[i] /= (GRID_X*GRID_Y);
            u_y_phys[i] /= (GRID_X*GRID_Y);
        }

        // --- Optional: Output velocity field for this time step ---
        if (step % 1000 == 0)
        {
            std::vector<double> max_vel;
            for (int i = 0; i < GRID_X; i++)
            {
                double x = i * LENGTH_X / GRID_X;
                for (int j = 0; j < GRID_Y; j++)
                {
                    double y = j * LENGTH_Y / GRID_Y;
                    int idx = i * GRID_Y + j;
                    double v = sqrt(pow(u_x_phys[idx],2) + pow(u_y_phys[idx],2));
                    max_vel.push_back(v);
                }                
            }
            std::cout << "Step: " << step << std::endl;
            double t = step * T_STEP;
            auto h = std::max_element(max_vel.begin(), max_vel.end());
            file_velocity_max << t << " " <<  *h << "\n";
            for (int i = 0; i < GRID_X; ++i) {
                double x = i * LENGTH_X / GRID_X;
                for (int j = 0; j < GRID_Y; ++j) {
                    double y = j * LENGTH_Y / GRID_Y;
                    int idx = i * GRID_Y + j;
                    file_velocity << x << " " << y << " " << sqrt(pow(u_x_phys[idx],2) + pow(u_y_phys[idx],2)) << " " << u_x_phys[idx] << " " <<  u_y_phys[idx]  << "\n";
                }
                file_velocity << "\n";
            }
            file_velocity << "\n" << std::flush;
        }

        // --- Step 4: Compute gradients of omega in Fourier space ---
        // Compute derivative fields:
        //   d(omega)/dx: multipLENGTH_Y omega_hat by i*kx.
        //   d(omega)/dy: multipLENGTH_Y omega_hat by i*ky.
        for (int i = 0; i < GRID_X; ++i) {
            for (int j = 0; j < GRID_Y/2+1; ++j) {
                int idx = i*(GRID_Y/2+1) + j;
                double reOmega = omega_hat[idx][0];
                double imOmega = omega_hat[idx][1];
                // d(omega)/dx_hat = i*kx*omega_hat => (-kx*im, kx*re)
                ux_hat[idx][0] = -kx[i]*imOmega;
                ux_hat[idx][1] =  kx[i]*reOmega;
                // d(omega)/dy_hat = i*ky*omega_hat => (-ky*im, ky*re)
                uy_hat[idx][0] = -ky[j]*imOmega;
                uy_hat[idx][1] =  ky[j]*reOmega;
            }
        }
        // Inverse FFT to recover gradients in physical space.
        fftw_plan plan_backward_gradx = fftw_plan_dft_c2r_2d(GRID_X, GRID_Y, ux_hat, omega, FFTW_MEASURE);
        fftw_plan plan_backward_grady = fftw_plan_dft_c2r_2d(GRID_X, GRID_Y, uy_hat, N_phys, FFTW_MEASURE);
        fftw_execute(plan_backward_gradx);
        fftw_execute(plan_backward_grady);
        // Store gradients in separate arrays for clarity.
        // We reuse omega array temporariLENGTH_Y for d(omega)/dx and N_phys for d(omega)/dy.
        std::vector<double> domega_dx(real_size), domega_dy(real_size);
        for (int i = 0; i < real_size; ++i) {
            domega_dx[i] = omega[i] / (GRID_X*GRID_Y);
            domega_dy[i] = N_phys[i] / (GRID_X*GRID_Y);
        }
        fftw_destroy_plan(plan_backward_gradx);
        fftw_destroy_plan(plan_backward_grady);

        // --- Step 5: Compute the nonlinear term in physical space ---
        // N_phys = u_x * d(omega)/dx + u_y * d(omega)/dy.
        for (int i = 0; i < real_size; ++i) {
            N_phys[i] = u_x_phys[i] * domega_dx[i] + u_y_phys[i] * domega_dy[i];
        }

        // --- Step 6: Transform the nonlinear term to Fourier space ---
        fftw_plan plan_forward_N = fftw_plan_dft_r2c_2d(GRID_X, GRID_Y, N_phys, N_hat, FFTW_MEASURE);
        fftw_execute(plan_forward_N);
        fftw_destroy_plan(plan_forward_N);

        // --- Step 7: Update omega_hat using an explicit Euler step ---
        // Equation in Fourier space:
        // d(omega_hat)/T_STEP = -N_hat - NU * |k|^2 * omega_hat.
        for (int i = 0; i < GRID_X; ++i) {
            for (int j = 0; j < GRID_Y/2+1; ++j) {
                int idx = i*(GRID_Y/2+1) + j;
                double k2 = kx[i]*kx[i] + ky[j]*ky[j];
                double rhs_re = -N_hat[idx][0] - NU * k2 * omega_hat[idx][0];
                double rhs_im = -N_hat[idx][1] - NU * k2 * omega_hat[idx][1];
                omega_hat[idx][0] += T_STEP * rhs_re;
                omega_hat[idx][1] += T_STEP * rhs_im;
            }
        }

        // --- Step 8: Inverse FFT to update omega in physical space ---
        fftw_execute(plan_backward_omega);
        for (int i = 0; i < real_size; ++i)
            omega[i] /= (GRID_X * GRID_Y);

        
        ens = 0;
        for (int i = 0; i < GRID_X; ++i) {
            double x = i * LENGTH_X / GRID_X;
            for (int j = 0; j < GRID_Y; ++j) {
                double y = j * LENGTH_Y / GRID_Y;
                int idx = i * GRID_Y + j;
                ens += sqrt(pow(u_x_phys[idx],2) + pow(u_y_phys[idx],2));
            }
        }
        ens *= (LENGTH_X/GRID_X) * (LENGTH_Y/GRID_Y);
        if (step % 1000 == 0 && step > 0)
        {
            file_ens << (step * T_STEP) << " " << (-0.5*(ens-ens_old)/(T_STEP*NU)) << "\n";
        }
        
    }

    // CleaNUp: Destroy FFTW plans and free allocated memory.
    fftw_destroy_plan(plan_forward_omega);
    fftw_destroy_plan(plan_backward_omega);
    fftw_destroy_plan(plan_backward_A);
    fftw_destroy_plan(plan_backward_ux);
    fftw_destroy_plan(plan_backward_uy);

    fftw_free(omega);
    fftw_free(N_phys);
    fftw_free(A_phys);
    fftw_free(omega_hat);
    fftw_free(A_hat);
    fftw_free(N_hat);
    fftw_free(u_x_phys);
    fftw_free(u_y_phys);
    fftw_free(ux_hat);
    fftw_free(uy_hat);
    file_velocity.close();
    file_velocity_max.close();
    file_ens.close();
    std::cout << "--- Normal Program Exit ---" << std::endl;

    return 0;
}
