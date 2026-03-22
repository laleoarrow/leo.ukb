(function () {
  const storedTheme = localStorage.getItem("theme");
  if (!storedTheme) {
    localStorage.setItem("theme", "dark");
    document.documentElement.setAttribute("data-bs-theme", "dark");
  }
})();
