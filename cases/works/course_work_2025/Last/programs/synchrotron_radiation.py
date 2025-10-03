# synchrotron_radiation.py
import numpy as np
import matplotlib.pyplot as plt
from Constants import C, M_SUN, PC, ERG, YR

class SynchrotronEmission:
    """Модель синхротронного излучения"""
    
    def __init__(self, p=2.2, epsilon_B=0.1, epsilon_e=0.1):
        self.p = p
        self.epsilon_B = epsilon_B
        self.epsilon_e = epsilon_e
        self.e = 4.803e-10
        self.m_e = 9.109e-28
        self.c = C
        self.sigma_T = 6.65e-25
    
    def magnetic_field(self, rho, v_shock):
        """Магнитное поле"""
        u_shock = 0.5 * rho * v_shock**2
        B = np.sqrt(8 * np.pi * self.epsilon_B * u_shock)
        return B
    
    def electron_lorentz_factor(self, v_shock):
        """Характерный лоренц-фактор электронов"""
        gamma_min = self.epsilon_e * (self.m_e * v_shock**2) / (2 * self.m_e * self.c**2)
        return gamma_min
    
    def critical_frequency(self, gamma, B):
        """Критическая частота синхротронного излучения"""
        nu_c = (3 * self.e * B * gamma**2) / (4 * np.pi * self.m_e * self.c)
        return nu_c
    
    def calculate_light_curve(self, times, radii, densities, velocities, nu=8.4e9, D=1e6*PC):
        """Расчет кривой блеска"""
        fluxes = []
        magnetic_fields = []
        lorentz_factors = []
        critical_frequencies = []
        
        for t, R, rho, v in zip(times, radii, densities, velocities):
            B = self.magnetic_field(rho, v)
            gamma_min = self.electron_lorentz_factor(v)
            nu_c = self.critical_frequency(gamma_min, B)
            
            # Упрощенный расчет потока
            F_nu = 1e-26 * (R/PC)**3 * (B/1e-4)**(self.p+1)/2 * (nu/1e9)**(-self.p/2) * (t/YR)**(-1)
            
            fluxes.append(F_nu)
            magnetic_fields.append(B)
            lorentz_factors.append(gamma_min)
            critical_frequencies.append(nu_c)
        
        return (np.array(fluxes), np.array(magnetic_fields), 
                np.array(lorentz_factors), np.array(critical_frequencies))

def plot_synchrotron_light_curves(synch, times, radii, densities, velocities, frequencies=[1.4e9, 8.4e9, 22e9]):
    """Построение кривых блеска на разных частотах"""
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Синхротронное излучение остатка сверхновой', fontsize=16, y=0.98)
    
    colors = ['blue', 'red', 'green']
    labels = [f'{nu/1e9:.1f} ГГц' for nu in frequencies]
    
    # Кривые блеска на разных частотах
    for idx, nu in enumerate(frequencies):
        fluxes, B_fields, gamma_min, nu_c = synch.calculate_light_curve(
            times, radii, densities, velocities, nu=nu)
        
        axes[0, 0].loglog(times/YR, fluxes, 
                         color=colors[idx], linewidth=2, 
                         label=labels[idx], marker='o', markersize=4)
    
    axes[0, 0].set_xlabel('Время (лет)')
    axes[0, 0].set_ylabel('Поток (эрг/с/см²/Гц)')
    axes[0, 0].set_title('Кривые блеска на разных частотах')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Эволюция магнитного поля
    fluxes, B_fields, gamma_min, nu_c = synch.calculate_light_curve(
        times, radii, densities, velocities)
    
    axes[0, 1].semilogy(times/YR, B_fields, 'b-', linewidth=2, marker='s', markersize=4)
    axes[0, 1].set_xlabel('Время (лет)')
    axes[0, 1].set_ylabel('Магнитное поле (Гс)')
    axes[0, 1].set_title('Эволюция магнитного поля')
    axes[0, 1].grid(True, alpha=0.3)
    
    # Эволюция лоренц-фактора электронов
    axes[1, 0].semilogy(times/YR, gamma_min, 'r-', linewidth=2, marker='^', markersize=4)
    axes[1, 0].set_xlabel('Время (лет)')
    axes[1, 0].set_ylabel('γ_min')
    axes[1, 0].set_title('Минимальный лоренц-фактор электронов')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Эволюция критической частоты
    axes[1, 1].semilogy(times/YR, nu_c/1e9, 'g-', linewidth=2, marker='d', markersize=4)
    axes[1, 1].set_xlabel('Время (лет)')
    axes[1, 1].set_ylabel('ν_c (ГГц)')
    axes[1, 1].set_title('Критическая частота синхротронного излучения')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('synchrotron_light_curves.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_frequency_spectrum(synch, times_years, radii, densities, velocities):
    """Построение спектров в разные моменты времени"""
    
    frequencies = np.logspace(8, 12, 100)  # от 100 МГц до 1000 ГГц
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    
    colors = ['blue', 'red', 'green']
    markers = ['o', 's', '^']
    
    for idx, t_years in enumerate(times_years):
        # Находим индекс ближайшего времени в массиве times
        t_target = t_years * YR
        t_idx = np.argmin(np.abs(times_years - t_years))
        
        # Проверяем, что индекс в пределах массива
        if t_idx >= len(radii):
            t_idx = len(radii) - 1
        
        fluxes_at_time = []
        for nu in frequencies:
            fluxes, _, _, _ = synch.calculate_light_curve(
                [t_target], [radii[t_idx]], [densities[t_idx]], 
                [velocities[t_idx]], nu=nu)
            fluxes_at_time.append(fluxes[0])
        
        ax.loglog(frequencies/1e9, fluxes_at_time, 
                 color=colors[idx], linewidth=2, 
                 label=f't = {t_years} лет', marker=markers[idx], markersize=4, alpha=0.7)
    
    ax.set_xlabel('Частота (ГГц)')
    ax.set_ylabel('Поток (эрг/с/см²/Гц)')
    ax.set_title('Спектры синхротронного излучения в разные моменты времени')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('synchrotron_spectra.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_parameter_dependence(synch_base):
    """Исследование зависимости от параметров ε_B и ε_e"""
    
    times = np.logspace(0, 3, 20) * YR
    radii = 2 * PC * (times/YR)**0.4
    densities = 1e-24 * np.ones_like(times)
    velocities = 1e9 * (times/YR)**(-0.6)
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Зависимость синхротронного излучения от параметров', fontsize=16)
    
    # Зависимость от epsilon_B
    epsilon_B_values = [0.01, 0.1, 0.3]
    colors = ['blue', 'red', 'green']
    
    for idx, eps_B in enumerate(epsilon_B_values):
        synch = SynchrotronEmission(epsilon_B=eps_B, epsilon_e=synch_base.epsilon_e)
        fluxes, B_fields, _, _ = synch.calculate_light_curve(times, radii, densities, velocities)
        
        axes[0, 0].loglog(times/YR, fluxes, 
                         color=colors[idx], linewidth=2, 
                         label=f'ε_B = {eps_B}')
        axes[0, 1].semilogy(times/YR, B_fields, 
                          color=colors[idx], linewidth=2, 
                          label=f'ε_B = {eps_B}')
    
    axes[0, 0].set_xlabel('Время (лет)')
    axes[0, 0].set_ylabel('Поток (эрг/с/см²/Гц)')
    axes[0, 0].set_title('Зависимость потока от ε_B (ν = 8.4 ГГц)')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    axes[0, 1].set_xlabel('Время (лет)')
    axes[0, 1].set_ylabel('Магнитное поле (Гс)')
    axes[0, 1].set_title('Зависимость магнитного поля от ε_B')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Зависимость от epsilon_e
    epsilon_e_values = [0.05, 0.1, 0.2]
    
    for idx, eps_e in enumerate(epsilon_e_values):
        synch = SynchrotronEmission(epsilon_B=synch_base.epsilon_B, epsilon_e=eps_e)
        fluxes, _, gamma_min, nu_c = synch.calculate_light_curve(times, radii, densities, velocities)
        
        axes[1, 0].loglog(times/YR, fluxes, 
                         color=colors[idx], linewidth=2, 
                         label=f'ε_e = {eps_e}')
        axes[1, 1].semilogy(times/YR, gamma_min, 
                          color=colors[idx], linewidth=2, 
                          label=f'ε_e = {eps_e}')
    
    axes[1, 0].set_xlabel('Время (лет)')
    axes[1, 0].set_ylabel('Поток (эрг/с/см²/Гц)')
    axes[1, 0].set_title('Зависимость потока от ε_e (ν = 8.4 ГГц)')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    axes[1, 1].set_xlabel('Время (лет)')
    axes[1, 1].set_ylabel('γ_min')
    axes[1, 1].set_title('Зависимость γ_min от ε_e')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('synchrotron_parameter_dependence.png', dpi=300, bbox_inches='tight')
    plt.show()

def generate_test_data():
    """Генерация тестовых данных для моделирования"""
    times = np.logspace(0, 3, 50) * YR  # от 1 до 1000 лет
    radii = 1.0 * PC * (times/YR)**0.4  # примерный закон расширения
    densities = 1e-24 * (times/YR)**(-0.5)  # примерная эволюция плотности
    velocities = 1e9 * (times/YR)**(-0.6)  # примерная эволюция скорости
    
    return times, radii, densities, velocities

def plot_frequency_spectrum_corrected(synch, times, radii, densities, velocities, selected_times_years=[10, 100, 1000]):
    """Исправленная версия построения спектров"""
    
    frequencies = np.logspace(8, 12, 100)  # от 100 МГц до 1000 ГГц
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    
    colors = ['blue', 'red', 'green']
    markers = ['o', 's', '^']
    
    for idx, t_years in enumerate(selected_times_years):
        # Находим индекс ближайшего времени в массиве times
        t_target = t_years * YR
        t_idx = np.argmin(np.abs(times/YR - t_years))
        
        # Проверяем, что индекс в пределах массива
        if t_idx >= len(radii):
            t_idx = len(radii) - 1
        
        print(f"Построение спектра для t = {t_years} лет (индекс {t_idx})")
        
        fluxes_at_time = []
        for nu in frequencies:
            fluxes, _, _, _ = synch.calculate_light_curve(
                [times[t_idx]], [radii[t_idx]], [densities[t_idx]], 
                [velocities[t_idx]], nu=nu)
            fluxes_at_time.append(fluxes[0])
        
        ax.loglog(frequencies/1e9, fluxes_at_time, 
                 color=colors[idx], linewidth=2, 
                 label=f't = {t_years} лет', marker=markers[idx], markersize=4, alpha=0.7)
    
    ax.set_xlabel('Частота (ГГц)')
    ax.set_ylabel('Поток (эрг/с/см²/Гц)')
    ax.set_title('Спектры синхротронного излучения в разные моменты времени')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('synchrotron_spectra.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    print("=" * 70)
    print("МОДЕЛИРОВАНИЕ СИНХРОТРОННОГО ИЗЛУЧЕНИЯ ОСТАТКОВ СВЕРХНОВЫХ")
    print("=" * 70)
    
    # Создаем объект модели синхротронного излучения
    synch = SynchrotronEmission(p=2.2, epsilon_B=0.1, epsilon_e=0.1)
    
    print("Параметры модели:")
    print(f"Степень распределения электронов: p = {synch.p}")
    print(f"Доля энергии в магнитном поле: ε_B = {synch.epsilon_B}")
    print(f"Доля энергии в электронах: ε_e = {synch.epsilon_e}")
    
    # Генерируем тестовые данные
    times, radii, densities, velocities = generate_test_data()
    
    # 1. Кривые блеска и эволюция параметров
    print("\n1. Построение кривых блеска и эволюции параметров...")
    plot_synchrotron_light_curves(synch, times, radii, densities, velocities)
    
    # 2. Спектры в разные моменты времени (исправленная версия)
    print("\n2. Построение спектров в разные моменты времени...")
    plot_frequency_spectrum_corrected(synch, times, radii, densities, velocities)
    
    # 3. Зависимость от параметров
    print("\n3. Исследование зависимости от параметров ε_B и ε_e...")
    plot_parameter_dependence(synch)
    
    # Расчет и вывод конкретных значений
    print("\n" + "=" * 80)
    print("ПАРАМЕТРЫ СИНХРОТРОННОГО ИЗЛУЧЕНИЯ:")
    print("=" * 80)
    print(f"{'Время (лет)':<12} {'B (мкГс)':<15} {'γ_min':<15} {'ν_c (ГГц)':<15} {'F_8.4ГГц (мЯ)':<15}")
    print("-" * 80)
    
    test_times_years = [1, 10, 100, 1000]
    for t_years in test_times_years:
        t_idx = np.argmin(np.abs(times/YR - t_years))
        
        fluxes, B_fields, gamma_min, nu_c = synch.calculate_light_curve(
            [times[t_idx]], [radii[t_idx]], [densities[t_idx]], [velocities[t_idx]])
        
        print(f"{t_years:<12} {B_fields[0]/1e-6:<15.2f} {gamma_min[0]:<15.2e} {nu_c[0]/1e9:<15.2f} {fluxes[0]*1e3:<15.2e}")
    
    print("\n" + "=" * 70)
    print("ГРАФИКИ СОХРАНЕНЫ:")
    print("- synchrotron_light_curves.png (кривые блеска и параметры)")
    print("- synchrotron_spectra.png (спектры в разные времена)")
    print("- synchrotron_parameter_dependence.png (зависимость от параметров)")
    print("=" * 70)
    