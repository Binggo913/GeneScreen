#!/usr/bin/env python3
"""Frontend assets embedded into the multi-query HTML report."""


LINKVIEW_SVG_HELPERS_JS = r"""
    const SVG_DEFINITION_TAGS = new Set([
      'defs',
      'clippath',
      'mask',
      'pattern',
      'symbol',
      'marker',
      'filter',
      'lineargradient',
      'radialgradient',
    ]);
    function isSvgDefinitionElement(el) {
      let node = el ? el.parentNode : null;
      const owner = el ? el.ownerSVGElement : null;
      while (node && node !== owner) {
        const tag = String(node.localName || node.nodeName || '').toLowerCase();
        if (SVG_DEFINITION_TAGS.has(tag)) return true;
        node = node.parentNode;
      }
      return false;
    }
    function isDisplaySuppressed(el) {
      const display = String(el.getAttribute('display') || '').toLowerCase();
      const visibility = String(el.getAttribute('visibility') || '').toLowerCase();
      const style = String(el.getAttribute('style') || '').toLowerCase().replace(/\s+/g, '');
      return display === 'none'
        || visibility === 'hidden'
        || style.includes('display:none')
        || style.includes('visibility:hidden');
    }
    function visibleChroRects(svg) {
      if (!svg) return [];
      return Array.from(svg.querySelectorAll('rect.chro'))
        .filter(rect => !isSvgDefinitionElement(rect))
        .filter(rect => !isDisplaySuppressed(rect))
        .sort((a, b) => Number(a.getAttribute('y') || 0) - Number(b.getAttribute('y') || 0));
    }
"""
