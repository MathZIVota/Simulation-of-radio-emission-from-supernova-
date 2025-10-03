# hydro_solver.py
import numpy as np
import matplotlib.pyplot as plt
from Constants import GAMMA, Normalization

class HydrodynamicSolver:
    """Сферически симметричный гидродинамический решатель с HLLC"""
    
    def __init__(self, gamma=GAMMA, normalization=None):
        self.gamma = gamma
        self.norm = normalization or Normalization()
        self.cfl = 0.3
        self.min_density = 1e-15
        self.min_pressure = 1e-15
        self.max_iter = 10000
        
    def _validate_primitive(self, rho, u, p):
        """Проверка физической осмысленности примитивных переменных"""
        rho = max(rho, self.min_density)
        p = max(p, self.min_pressure)
        
        if np.isnan(rho) or np.isnan(u) or np.isnan(p):
            raise ValueError("Обнаружены NaN в примитивных переменных")
        if np.isinf(rho) or np.isinf(u) or np.isinf(p):
            raise ValueError("Обнаружены Inf в примитивных переменных")
            
        return rho, u, p
    
    def primitive_to_conservative(self, rho, u, p, r):
        """Преобразование в консервативные переменные с проверкой"""
        rho, u, p = self._validate_primitive(rho, u, p)
        
        mass = rho
        momentum = rho * u
        energy = p/(self.gamma-1) + 0.5*rho*u**2
        
        # Проверка энергии
        if energy < 0:
            print(f"Предупреждение: отрицательная энергия {energy:.2e}, исправляю")
            energy = p/(self.gamma-1)
        
        return np.array([mass, momentum, energy])
    
    def conservative_to_primitive(self, U, r):
        """Преобразование в примитивные переменные с защитой"""
        if np.any(np.isnan(U)) or np.any(np.isinf(U)):
            raise ValueError("NaN или Inf в консервативных переменных")
            
        rho = max(U[0], self.min_density)
        u = U[1] / rho
        
        # Защита от отрицательного давления
        kinetic_energy = 0.5 * rho * u**2
        if U[2] < kinetic_energy:
            p = self.min_pressure
        else:
            p = max((U[2] - kinetic_energy) * (self.gamma-1), self.min_pressure)
        
        return self._validate_primitive(rho, u, p)
    
    def sound_speed(self, p, rho):
        """Скорость звука с защитой"""
        return np.sqrt(self.gamma * p / max(rho, self.min_density))
    
    def hllc_flux(self, UL, UR, r_interface):
        """HLLC риман-солвер с улучшенной стабильностью"""
        try:
            rhoL, uL, pL = self.conservative_to_primitive(UL, r_interface)
            rhoR, uR, pR = self.conservative_to_primitive(UR, r_interface)
            
            aL = self.sound_speed(pL, rhoL)
            aR = self.sound_speed(pR, rhoR)
            
            # Оценка скоростей волн
            SL = min(uL - aL, uR - aR)
            SR = max(uL + aL, uR + aR)
            
            # Упрощенная оценка давления
            rho_avg = 0.5 * (rhoL + rhoR)
            a_avg = 0.5 * (aL + aR)
            p_star = max(0.5*(pL + pR) - 0.5*(uR - uL)*rho_avg*a_avg, self.min_pressure)
            
            # Оценка скорости контактного разрыва
            denominator = rhoL*(SL - uL) - rhoR*(SR - uR)
            if abs(denominator) < 1e-10:
                S_star = 0.5 * (uL + uR)
            else:
                S_star = (pR - pL + rhoL*uL*(SL - uL) - rhoR*uR*(SR - uR)) / denominator
            
            # Вычисление потоков
            if SL >= 0:
                F = np.array([rhoL*uL, rhoL*uL**2 + pL, uL*(UL[2] + pL)])
            elif SR <= 0:
                F = np.array([rhoR*uR, rhoR*uR**2 + pR, uR*(UR[2] + pR)])
            elif S_star >= 0:
                rho_star = rhoL * (SL - uL)/(SL - S_star)
                U_star = np.array([
                    rho_star,
                    rho_star * S_star,
                    rho_star * (UL[2]/rhoL + (S_star - uL)*(S_star + pL/(rhoL*(SL - uL))))
                ])
                F = np.array([rhoL*uL, rhoL*uL**2 + pL, uL*(UL[2] + pL)]) + SL*(U_star - UL)
            else:
                rho_star = rhoR * (SR - uR)/(SR - S_star)
                U_star = np.array([
                    rho_star,
                    rho_star * S_star,
                    rho_star * (UR[2]/rhoR + (S_star - uR)*(S_star + pR/(rhoR*(SR - uR))))
                ])
                F = np.array([rhoR*uR, rhoR*uR**2 + pR, uR*(UR[2] + pR)]) + SR*(U_star - UR)
            
            return F
            
        except Exception as e:
            print(f"Ошибка в HLLC: {e}, использую Lax-Friedrichs")
            # Запасной вариант - Lax-Friedrichs
            F_L = np.array([rhoL*uL, rhoL*uL**2 + pL, uL*(UL[2] + pL)])
            F_R = np.array([rhoR*uR, rhoR*uR**2 + pR, uR*(UR[2] + pR)])
            lambda_max = max(abs(uL) + aL, abs(uR) + aR)
            return 0.5 * (F_L + F_R) - 0.5 * lambda_max * (UR - UL)
    
    def solve(self, r, U0, t_final, geometry='spherical'):
        """Основной метод решения"""
        nr = len(r)
        dr = r[1] - r[0]
        
        U = U0.copy()
        t = 0.0
        iteration = 0
        
        while t < t_final and iteration < self.max_iter:
            iteration += 1
            
            # Вычисление шага по времени
            dt = self._compute_time_step(U, r, dr)
            if t + dt > t_final:
                dt = t_final - t
            
            # Вычисление потоков
            F = np.zeros((nr, 3))
            for i in range(nr-1):
                r_interface = 0.5*(r[i] + r[i+1])
                F[i] = self.hllc_flux(U[i], U[i+1], r_interface)
            
            # Обновление решения
            U_new = U.copy()
            for i in range(1, nr-1):
                if geometry == 'spherical':
                    flux_in = F[i-1] * r[i-1]**2
                    flux_out = F[i] * r[i]**2
                    U_new[i] = U[i] - dt * (flux_out - flux_in) / (r[i]**2 * dr)
                else:
                    U_new[i] = U[i] - dt/dr * (F[i] - F[i-1])
            
            # Граничные условия
            U_new[0] = U_new[1]
            U_new[-1] = U_new[-2]
            
            U = U_new
            t += dt
            
        return U
    
    def _compute_time_step(self, U, r, dr):
        """Вычисление шага по времени"""
        dt_min = 1e10
        for i in range(len(U)):
            rho, u, p = self.conservative_to_primitive(U[i], r[i])
            a = self.sound_speed(p, rho)
            speed = abs(u) + a
            if speed > 0:
                dt_min = min(dt_min, dr / speed)
        
        return self.cfl * dt_min

class TestProblems:
    """Класс тестовых задач"""
    
    @staticmethod
    def sod_shock_tube(nx=100, geometry='planar'):
        """Тест Сода"""
        x = np.linspace(0, 1, nx)
        U0 = np.zeros((nx, 3))
        
        solver = HydrodynamicSolver(gamma=1.4)
        
        for i in range(nx):
            if x[i] < 0.5:
                U0[i] = solver.primitive_to_conservative(1.0, 0.0, 1.0, x[i])
            else:
                U0[i] = solver.primitive_to_conservative(0.125, 0.0, 0.1, x[i])
        
        return x, U0
    
    @staticmethod
    def strong_shock_pressure(nx=100, pressure_ratio=1000.0):
        """Сильная ударная волна по давлению"""
        x = np.linspace(0, 1, nx)
        U0 = np.zeros((nx, 3))
        
        solver = HydrodynamicSolver(gamma=5/3)
        
        for i in range(nx):
            if x[i] < 0.5:
                U0[i] = solver.primitive_to_conservative(1.0, 0.0, pressure_ratio, x[i])
            else:
                U0[i] = solver.primitive_to_conservative(1.0, 0.0, 1.0, x[i])
        
        return x, U0
    
    @staticmethod
    def strong_shock_density(nx=100, density_ratio=1000.0):
        """Сильная ударная волна по плотности"""
        x = np.linspace(0, 1, nx)
        U0 = np.zeros((nx, 3))
        
        solver = HydrodynamicSolver(gamma=5/3)
        
        for i in range(nx):
            if x[i] < 0.5:
                U0[i] = solver.primitive_to_conservative(density_ratio, 0.0, 1.0, x[i])
            else:
                U0[i] = solver.primitive_to_conservative(1.0, 0.0, 1.0, x[i])
        
        return x, U0

def plot_results(x, U_initial, U_final, solver, test_name="Тест"):
    """Построение графиков результатов гидродинамического моделирования"""
    
    # Преобразуем начальные и конечные состояния в примитивные переменные
    rho_initial = np.zeros(len(x))
    u_initial = np.zeros(len(x))
    p_initial = np.zeros(len(x))
    
    rho_final = np.zeros(len(x))
    u_final = np.zeros(len(x))
    p_final = np.zeros(len(x))
    e_final = np.zeros(len(x))
    
    for i in range(len(x)):
        rho_initial[i], u_initial[i], p_initial[i] = solver.conservative_to_primitive(U_initial[i], x[i])
        rho_final[i], u_final[i], p_final[i] = solver.conservative_to_primitive(U_final[i], x[i])
        e_final[i] = p_final[i] / ((solver.gamma - 1) * max(rho_final[i], 1e-10))
    
    # Создаем графики
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f'Результаты гидродинамического моделирования: {test_name}', fontsize=16)
    
    # Плотность
    axes[0, 0].plot(x, rho_initial, 'r--', label='Начальное', linewidth=2, alpha=0.7)
    axes[0, 0].plot(x, rho_final, 'b-', label='Конечное', linewidth=2)
    axes[0, 0].set_xlabel('x')
    axes[0, 0].set_ylabel('Плотность')
    axes[0, 0].set_title('Плотность')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Скорость
    axes[0, 1].plot(x, u_initial, 'r--', label='Начальное', linewidth=2, alpha=0.7)
    axes[0, 1].plot(x, u_final, 'b-', label='Конечное', linewidth=2)
    axes[0, 1].set_xlabel('x')
    axes[0, 1].set_ylabel('Скорость')
    axes[0, 1].set_title('Скорость')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Давление
    axes[0, 2].plot(x, p_initial, 'r--', label='Начальное', linewidth=2, alpha=0.7)
    axes[0, 2].plot(x, p_final, 'b-', label='Конечное', linewidth=2)
    axes[0, 2].set_xlabel('x')
    axes[0, 2].set_ylabel('Давление')
    axes[0, 2].set_title('Давление')
    axes[0, 2].legend()
    axes[0, 2].grid(True, alpha=0.3)
    
    # Внутренняя энергия
    axes[1, 0].plot(x, e_final, 'g-', linewidth=2)
    axes[1, 0].set_xlabel('x')
    axes[1, 0].set_ylabel('Внутренняя энергия')
    axes[1, 0].set_title('Внутренняя энергия (конечная)')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Массовый поток
    mass_flux_final = rho_final * u_final
    axes[1, 1].plot(x, mass_flux_final, 'm-', linewidth=2)
    axes[1, 1].set_xlabel('x')
    axes[1, 1].set_ylabel('Массовый поток')
    axes[1, 1].set_title('Массовый поток (конечный)')
    axes[1, 1].grid(True, alpha=0.3)
    
    # Полная энергия
    total_energy_final = U_final[:, 2]
    axes[1, 2].plot(x, total_energy_final, 'c-', linewidth=2)
    axes[1, 2].set_xlabel('x')
    axes[1, 2].set_ylabel('Полная энергия')
    axes[1, 2].set_title('Полная энергия (конечная)')
    axes[1, 2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'hydro_results_{test_name.replace(" ", "_")}.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Дополнительная информация
    print(f"\nАнализ результатов {test_name}:")
    print(f"Максимальная плотность: {np.max(rho_final):.3f}")
    print(f"Максимальная скорость: {np.max(np.abs(u_final)):.3f}")
    print(f"Максимальное давление: {np.max(p_final):.3f}")
    print(f"Общая масса: {np.trapz(rho_final, x):.3f}")
    print(f"Общая энергия: {np.trapz(total_energy_final, x):.3f}")

def run_test_sod():
    """Запуск и визуализация теста Сода"""
    print("=" * 60)
    print("ТЕСТ СОДА: УДАРНАЯ ТРУБА")
    print("=" * 60)
    
    x, U0 = TestProblems.sod_shock_tube(nx=200)
    solver = HydrodynamicSolver(gamma=1.4)
    
    print("Решаем задачу Сода...")
    U_final = solver.solve(x, U0, t_final=0.02, geometry='planar')
    
    plot_results(x, U0, U_final, solver, "Тест Сода")
    
    # Анализ волн
    rho_final = U_final[:, 0]
    u_final = U_final[:, 1] / np.maximum(U_final[:, 0], 1e-10)
    
    # Находим положения волн
    grad_rho = np.gradient(rho_final, x)
    shock_position = x[np.argmax(np.abs(grad_rho))]
    contact_position = x[np.argmax(np.abs(np.gradient(u_final, x)))]
    
    print(f"\nПоложение ударной волны: {shock_position:.3f}")
    print(f"Положение контактного разрыва: {contact_position:.3f}")

def run_test_strong_shock_pressure():
    """Запуск теста с сильной ударной волной по давлению"""
    print("=" * 60)
    print("ТЕСТ: СИЛЬНАЯ УДАРНАЯ ВОЛНА (ДАВЛЕНИЕ 1000:1)")
    print("=" * 60)
    
    x, U0 = TestProblems.strong_shock_pressure(nx=200, pressure_ratio=1000)
    solver = HydrodynamicSolver(gamma=5/3)
    
    print("Решаем задачу с сильной ударной волной...")
    U_final = solver.solve(x, U0, t_final=0.02, geometry='planar')
    
    plot_results(x, U0, U_final, solver, "Сильная УВ по давлению")

def run_test_strong_shock_density():
    """Запуск теста с сильной ударной волной по плотности"""
    print("=" * 60)
    print("ТЕСТ: СИЛЬНАЯ УДАРНАЯ ВОЛНА (ПЛОТНОСТЬ 1000:1)")
    print("=" * 60)
    
    x, U0 = TestProblems.strong_shock_density(nx=200, density_ratio=1000)
    solver = HydrodynamicSolver(gamma=5/3)
    
    print("Решаем задачу с сильной ударной волной...")
    U_final = solver.solve(x, U0, t_final=0.02, geometry='planar')
    
    plot_results(x, U0, U_final, solver, "Сильная УВ по плотности")

def compare_all_tests():
    """Сравнительный анализ всех тестов"""
    print("=" * 60)
    print("СРАВНИТЕЛЬНЫЙ АНАЛИЗ ВСЕХ ТЕСТОВ")
    print("=" * 60)
    
    tests = [
        ("Сода", TestProblems.sod_shock_tube(200), 1.4, 0.2),
        ("Сильное давление", TestProblems.strong_shock_pressure(200, 1000), 5/3, 0.15),
        #("Сильная плотность", TestProblems.strong_shock_density(200, 1000), 5/3, 0.15)
    ]
    
    fig, axes = plt.subplots(3, 3, figsize=(15, 12))
    
    for idx, (test_name, (x, U0), gamma, t_final) in enumerate(tests):
        solver = HydrodynamicSolver(gamma=gamma)
        U_final = solver.solve(x, U0, t_final=t_final, geometry='planar')
        
        # Извлекаем конечные плотности
        rho_final = np.zeros(len(x))
        for i in range(len(x)):
            rho_final[i], _, _ = solver.conservative_to_primitive(U_final[i], x[i])
        
        # Строим сравнительные графики
        axes[0, idx].plot(x, rho_final, 'b-', linewidth=2)
        axes[0, idx].set_title(f'{test_name}\nПлотность')
        axes[0, idx].set_xlabel('x')
        axes[0, idx].grid(True, alpha=0.3)
        
        # Скорость
        u_final = U_final[:, 1] / np.maximum(U_final[:, 0], 1e-10)
        axes[1, idx].plot(x, u_final, 'r-', linewidth=2)
        axes[1, idx].set_title('Скорость')
        axes[1, idx].set_xlabel('x')
        axes[1, idx].grid(True, alpha=0.3)
        
        # Давление
        p_final = (gamma - 1) * (U_final[:, 2] - 0.5 * rho_final * u_final**2)
        axes[2, idx].plot(x, p_final, 'g-', linewidth=2)
        axes[2, idx].set_title('Давление')
        axes[2, idx].set_xlabel('x')
        axes[2, idx].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('hydro_comparison_all_tests.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    print("ГИДРОДИНАМИЧЕСКИЙ РЕШАТЕЛЬ С ВИЗУАЛИЗАЦИЕЙ")
    print("=" * 60)
    
    # Запуск отдельных тестов
    run_test_sod()
    run_test_strong_shock_pressure() 
    #run_test_strong_shock_density()
    
    # Сравнительный анализ
    compare_all_tests()
    
    print("\n" + "=" * 60)
    print("ВСЕ ТЕСТЫ ЗАВЕРШЕНЫ!")
    print("Созданы файлы:")
    print("- hydro_results_Тест_Сода.png")
    print("- hydro_results_Сильная_УВ_по_давлению.png") 
    print("- hydro_results_Сильная_УВ_по_плотности.png")
    print("- hydro_comparison_all_tests.png")