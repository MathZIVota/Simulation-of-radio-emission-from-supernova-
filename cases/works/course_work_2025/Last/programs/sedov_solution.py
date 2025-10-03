# sedov_solution.py
import numpy as np
import matplotlib.pyplot as plt
from Constants import GAMMA, M_SUN, PC, YR, ERG, KB, MP, MU

class SedovTaylorSolution:
    """Аналитическое решение Седова-Тейлора"""
    
    def __init__(self, E0=1e51, rho0=1e-24, gamma=GAMMA):
        self.E0 = E0
        self.rho0 = rho0
        self.gamma = gamma
        self.alpha = self._compute_alpha()
    
    def _compute_alpha(self):
        """Вычисление константы alpha"""
        if abs(self.gamma - 5/3) < 1e-3:
            return 1.15
        elif abs(self.gamma - 7/5) < 1e-3:
            return 1.03
        else:
            return 1.15
    
    def shock_radius(self, t):
        """Радиус ударной волны"""
        return self.alpha * (self.E0 / self.rho0)**(1/5) * t**(2/5)
    
    def shock_velocity(self, t):
        """Скорость ударной волны"""
        return (2/5) * self.alpha * (self.E0 / self.rho0)**(1/5) * t**(-3/5)
    
    def density_profile(self, r, t):
        """Плотность"""
        R = self.shock_radius(t)
        if r > R:
            return self.rho0
        xi = r / R
        if xi < 0.5:
            return 4 * self.rho0
        elif xi < 0.8:
            return 3 * self.rho0
        else:
            return 6 * self.rho0
    
    def pressure_profile(self, r, t):
        """Давление"""
        R = self.shock_radius(t)
        if r > R:
            return 1e-10
        v_s = self.shock_velocity(t)
        p2 = 2/(self.gamma + 1) * self.rho0 * v_s**2
        xi = r / R
        return p2 * (1 - 0.5*xi**2)
    
    def velocity_profile(self, r, t):
        """Скорость"""
        R = self.shock_radius(t)
        if r > R:
            return 0.0
        v_s = self.shock_velocity(t)
        xi = r / R
        return 2/(self.gamma + 1) * v_s * xi
    
    def temperature_profile(self, r, t):
        """Температура (K)"""
        rho = self.density_profile(r, t)
        p = self.pressure_profile(r, t)
        return p * MU * MP / (rho * KB)

def plot_sedov_profiles_comparison(sedov, times_years=[10, 100, 1000]):
    """Построение профилей параметров для нескольких моментов времени на одних графиках"""
    
    # Создаем фигуру с subplots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Профили параметров остатка сверхновой по решению Седова-Тейлора', 
                 fontsize=16, y=0.98)
    
    colors = ['blue', 'red', 'green']
    line_styles = ['-', '--', '-.']
    labels = [f'{t} лет' for t in times_years]
    
    for idx, t_years in enumerate(times_years):
        t = t_years * YR
        R = sedov.shock_radius(t)
        
        # Создаем радиальную сетку
        r = np.linspace(0.01 * R, 1.5 * R, 300)
        
        # Вычисляем профили
        density = np.array([sedov.density_profile(ri, t) for ri in r])
        pressure = np.array([sedov.pressure_profile(ri, t) for ri in r])
        velocity = np.array([sedov.velocity_profile(ri, t) for ri in r])
        temperature = np.array([sedov.temperature_profile(ri, t) for ri in r])
        
        # Нормировочные величины
        p_shock = sedov.pressure_profile(0.99*R, t)
        v_shock = sedov.shock_velocity(t)
        
        # Плотность (нормированная)
        axes[0, 0].plot(r/PC, density/sedov.rho0, 
                       color=colors[idx], linestyle=line_styles[idx], 
                       linewidth=2, label=labels[idx])
        axes[0, 0].axvline(R/PC, color=colors[idx], linestyle=':', alpha=0.7)
        
        # Давление (нормированное)
        axes[0, 1].plot(r/PC, pressure/p_shock, 
                       color=colors[idx], linestyle=line_styles[idx], 
                       linewidth=2, label=labels[idx])
        axes[0, 1].axvline(R/PC, color=colors[idx], linestyle=':', alpha=0.7)
        
        # Скорость (нормированная)
        axes[1, 0].plot(r/PC, velocity/v_shock, 
                       color=colors[idx], linestyle=line_styles[idx], 
                       linewidth=2, label=labels[idx])
        axes[1, 0].axvline(R/PC, color=colors[idx], linestyle=':', alpha=0.7)
        
        # Температура
        axes[1, 1].semilogy(r/PC, temperature, 
                           color=colors[idx], linestyle=line_styles[idx],
                           linewidth=2, label=labels[idx])
        axes[1, 1].axvline(R/PC, color=colors[idx], linestyle=':', alpha=0.7)
    
    # Настройка подграфиков
    axes[0, 0].set_xlabel('Радиус (пк)')
    axes[0, 0].set_ylabel('ρ / ρ₀')
    axes[0, 0].set_title('Нормированная плотность')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    axes[0, 1].set_xlabel('Радиус (пк)')
    axes[0, 1].set_ylabel('P / P_shock')
    axes[0, 1].set_title('Нормированное давление')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    axes[1, 0].set_xlabel('Радиус (пк)')
    axes[1, 0].set_ylabel('v / v_shock')
    axes[1, 0].set_title('Нормированная скорость')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    axes[1, 1].set_xlabel('Радиус (пк)')
    axes[1, 1].set_ylabel('Температура (K)')
    axes[1, 1].set_title('Температура')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('sedov_profiles_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Вывод информации о параметрах
    print("\n" + "=" * 80)
    print("ПАРАМЕТРЫ УДАРНОЙ ВОЛНЫ В РАЗНЫЕ МОМЕНТЫ ВРЕМЕНИ:")
    print("=" * 80)
    print(f"{'Время (лет)':<12} {'Радиус (пк)':<15} {'Скорость (тыс. км/с)':<20} {'Темп. за фронтом (K)':<20}")
    print("-" * 80)
    
    for t_years in times_years:
        t = t_years * YR
        R = sedov.shock_radius(t)
        v = sedov.shock_velocity(t)
        T_shock = sedov.temperature_profile(0.9*R, t)
        
        print(f"{t_years:<12} {R/PC:<15.2f} {v/1e8:<20.2f} {T_shock:<20.2e}")

def plot_sedov_evolution(sedov):
    """Построение эволюции параметров со временем"""
    times = np.logspace(0, 3, 50) * YR  # от 1 до 1000 лет
    
    radii = []
    velocities = []
    
    for t in times:
        R = sedov.shock_radius(t)
        v = sedov.shock_velocity(t)
        radii.append(R)
        velocities.append(v)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Эволюция остатка сверхновой по решению Седова-Тейлора', fontsize=16)
    
    # Радиус ударной волны
    axes[0, 0].loglog(times/YR, np.array(radii)/PC, 'b-', linewidth=2)
    axes[0, 0].set_xlabel('Время (лет)')
    axes[0, 0].set_ylabel('Радиус (пк)')
    axes[0, 0].set_title('Радиус ударной волны')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].annotate(f'R ∝ t$^{{{2/5:.2f}}}$', xy=(0.7, 0.1), xycoords='axes fraction', 
                       fontsize=12, bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
    
    # Скорость ударной волны
    axes[0, 1].loglog(times/YR, np.array(velocities)/1e8, 'r-', linewidth=2)
    axes[0, 1].set_xlabel('Время (лет)')
    axes[0, 1].set_ylabel('Скорость (1000 км/с)')
    axes[0, 1].set_title('Скорость ударной волны')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].annotate(f'v ∝ t$^{{{-3/5:.2f}}}$', xy=(0.7, 0.1), xycoords='axes fraction', 
                       fontsize=12, bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
    
    # Кинетическая энергия
    kinetic_energy = 0.5 * sedov.rho0 * (4/3 * np.pi * np.array(radii)**3) * np.array(velocities)**2
    axes[1, 0].semilogx(times/YR, kinetic_energy/sedov.E0, 'g-', linewidth=2)
    axes[1, 0].set_xlabel('Время (лет)')
    axes[1, 0].set_ylabel('E_kin / E_0')
    axes[1, 0].set_title('Доля кинетической энергии')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Число Маха
    T_ism = 1e4  # K
    c_s_ism = np.sqrt(GAMMA * KB * T_ism / (MU * MP))
    mach_numbers = np.array(velocities) / c_s_ism
    axes[1, 1].loglog(times/YR, mach_numbers, 'm-', linewidth=2)
    axes[1, 1].set_xlabel('Время (лет)')
    axes[1, 1].set_ylabel('Число Маха')
    axes[1, 1].set_title('Число Маха ударной волны')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('sedov_evolution.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    return times, radii, velocities

def plot_energy_dependence(sedov_base):
    """Исследование зависимости от энергии взрыва"""
    energies = [1e50, 1e51, 1e52]  # эрг
    colors = ['blue', 'red', 'green']
    labels = [f'E = {e:.0e} эрг' for e in energies]
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle('Зависимость эволюции от энергии взрыва', fontsize=14)
    
    times = np.logspace(0, 3, 50) * YR
    
    for idx, E in enumerate(energies):
        sedov = SedovTaylorSolution(E0=E, rho0=sedov_base.rho0, gamma=sedov_base.gamma)
        
        radii = [sedov.shock_radius(t) for t in times]
        velocities = [sedov.shock_velocity(t) for t in times]
        
        axes[0].loglog(times/YR, np.array(radii)/PC, 
                      color=colors[idx], linewidth=2, label=labels[idx])
        axes[1].loglog(times/YR, np.array(velocities)/1e8, 
                      color=colors[idx], linewidth=2, label=labels[idx])
    
    axes[0].set_xlabel('Время (лет)')
    axes[0].set_ylabel('Радиус (пк)')
    axes[0].set_title('Радиус ударной волны')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    axes[1].set_xlabel('Время (лет)')
    axes[1].set_ylabel('Скорость (1000 км/с)')
    axes[1].set_title('Скорость ударной волны')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('sedov_energy_dependence.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_density_dependence(sedov_base):
    """Исследование зависимости от плотности среды"""
    densities = [1e-25, 1e-24, 1e-23]  # г/см³
    colors = ['blue', 'red', 'green']
    labels = [f'ρ = {rho:.0e} г/см³' for rho in densities]
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle('Зависимость эволюции от плотности среды', fontsize=14)
    
    times = np.logspace(0, 3, 50) * YR
    
    for idx, rho in enumerate(densities):
        sedov = SedovTaylorSolution(E0=sedov_base.E0, rho0=rho, gamma=sedov_base.gamma)
        
        radii = [sedov.shock_radius(t) for t in times]
        velocities = [sedov.shock_velocity(t) for t in times]
        
        axes[0].loglog(times/YR, np.array(radii)/PC, 
                      color=colors[idx], linewidth=2, label=labels[idx])
        axes[1].loglog(times/YR, np.array(velocities)/1e8, 
                      color=colors[idx], linewidth=2, label=labels[idx])
    
    axes[0].set_xlabel('Время (лет)')
    axes[0].set_ylabel('Радиус (пк)')
    axes[0].set_title('Радиус ударной волны')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    axes[1].set_xlabel('Время (лет)')
    axes[1].set_ylabel('Скорость (1000 км/с)')
    axes[1].set_title('Скорость ударной волны')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('sedov_density_dependence.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    print("=" * 70)
    print("РЕШЕНИЕ СЕДОВА-ТЕЙЛОРА ДЛЯ ОСТАТКОВ СВЕРХНОВЫХ")
    print("=" * 70)
    
    # Создаем объект решения Седова
    sedov = SedovTaylorSolution(E0=1e51, rho0=1e-24, gamma=5/3)
    
    print("Параметры модели:")
    print(f"Энергия взрыва: {sedov.E0:.1e} эрг")
    print(f"Плотность среды: {sedov.rho0:.1e} г/см³")
    print(f"Показатель адиабаты: γ = {sedov.gamma:.2f}")
    print(f"Константа Седова: α = {sedov.alpha:.2f}")
    
    # 1. Сравнение профилей на одних графиках
    print("\n1. Построение сравнения профилей для разных времен...")
    plot_sedov_profiles_comparison(sedov, times_years=[10, 100, 1000])
    
    # 2. Эволюция со временем
    print("\n2. Построение эволюции параметров...")
    times, radii, velocities = plot_sedov_evolution(sedov)
    
    # 3. Зависимость от энергии
    print("\n3. Исследование зависимости от энергии взрыва...")
    plot_energy_dependence(sedov)
    
    # 4. Зависимость от плотности
    print("\n4. Исследование зависимости от плотности среды...")
    plot_density_dependence(sedov)
    
    print("\n" + "=" * 70)
    print("ГРАФИКИ СОХРАНЕНЫ:")
    print("- sedov_profiles_comparison.png (все профили в одном окне)")
    print("- sedov_evolution.png")
    print("- sedov_energy_dependence.png")
    print("- sedov_density_dependence.png")
    print("=" * 70)